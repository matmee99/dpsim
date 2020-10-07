/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/DP/DP_Ph1_ControlledVoltageSource.h>

using namespace CPS;

DP::Ph1::ControlledVoltageSource::ControlledVoltageSource(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Complex>(uid, name, logLevel) {
	setVirtualNodeNumber(1);
	setTerminalNumber(2);
	mIntfVoltage = MatrixComp::Zero(1, 1);
	mIntfCurrent = MatrixComp::Zero(1, 1);
}

void DP::Ph1::ControlledVoltageSource::setParameters(MatrixComp voltageRefABC) {
	mIntfVoltage = voltageRefABC;
	mParametersSet = true;
}

SimPowerComp<Complex>::Ptr DP::Ph1::ControlledVoltageSource::clone(String name) {
	auto copy = ControlledVoltageSource::make(name, mLogLevel);
	copy->setParameters(attribute<MatrixComp>("v_intf")->get());
	return copy;
}

void DP::Ph1::ControlledVoltageSource::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {

	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void DP::Ph1::ControlledVoltageSource::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		if (terminalNotGrounded(0)) {
			Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0), Complex(-1, 0), mNumFreqs, freq);
			Math::setMatrixElement(systemMatrix, matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex(), Complex(-1, 0), mNumFreqs, freq);
		}
		if (terminalNotGrounded(1)) {
			Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1), Complex(1, 0), mNumFreqs, freq);
			Math::setMatrixElement(systemMatrix, matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex(), Complex(1, 0), mNumFreqs, freq);
		}

		mSLog->info("-- Stamp frequency {:d} ---", freq);
		if (terminalNotGrounded(0)) {
			mSLog->info("Add {:f} to system at ({:d},{:d})", -1., matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
			mSLog->info("Add {:f} to system at ({:d},{:d})", -1., mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
		}
		if (terminalNotGrounded(1)) {
			mSLog->info("Add {:f} to system at ({:d},{:d})", 1., mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
			mSLog->info("Add {:f} to system at ({:d},{:d})", 1., matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
		}
	}
}

void DP::Ph1::ControlledVoltageSource::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(), mIntfVoltage(0, 0), mNumFreqs);
	SPDLOG_LOGGER_DEBUG(mSLog, "Add {:s} to source vector at {:d}",
		Logger::complexToString(mIntfVoltage(0, 0)), mVirtualNodes[0]->matrixNodeIndex());
}

void DP::Ph1::ControlledVoltageSource::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	attributeDependencies.push_back(this->attribute("v_intf"));
	modifiedAttributes.push_back(this->attribute("right_vector"));
}

void DP::Ph1::ControlledVoltageSource::mnaPreStep(Real time, Int timeStepCount) {
	this->mnaApplyRightSideVectorStamp(this->mRightVector);
}

void DP::Ph1::ControlledVoltageSource::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(this->attribute("i_intf"));
}

void DP::Ph1::ControlledVoltageSource::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	this->mnaUpdateCurrent(*leftVector);
}

void DP::Ph1::ControlledVoltageSource::mnaUpdateCurrent(const Matrix& leftVector) {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		mIntfCurrent(0, freq) = Math::complexFromVectorElement(leftVector, mVirtualNodes[0]->matrixNodeIndex(), mNumFreqs, freq);
	}
}
