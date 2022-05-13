/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/DP/DP_Ph1_SynchronGeneratorVBR.h>

using namespace CPS;

DP::Ph1::SynchronGeneratorVBR::SynchronGeneratorVBR
    (String uid, String name, Logger::Level logLevel)
	: Base::ReducedOrderSynchronGenerator<Complex>(uid, name, logLevel) {

	setVirtualNodeNumber(2);
	setTerminalNumber(1);
	
	// model variables
	mIntfVoltage = MatrixComp::Zero(1, 1);
	mIntfCurrent = MatrixComp::Zero(1, 1);

    // Register attributes
    addAttribute<Complex>("Evbr", &mEvbr, Flags::read);

	// initialize conductance Matrix
    mConductanceMatrix = Matrix::Zero(2,2);

	//
	mShiftVector = Matrix::Zero(3,1);
	mShiftVector << Complex(1., 0), SHIFT_TO_PHASE_B, SHIFT_TO_PHASE_C;
	mShiftVectorConj = Matrix::Zero(3,1);
	mShiftVectorConj << Complex(1., 0), std::conj(SHIFT_TO_PHASE_B), std::conj(SHIFT_TO_PHASE_C);
}

DP::Ph1::SynchronGeneratorVBR::SynchronGeneratorVBR
	(String name, Logger::Level logLevel)
	: SynchronGeneratorVBR(name, name, logLevel) {
}

DP::Ph1::SynchronGeneratorVBR::~SynchronGeneratorVBR() {
}

void DP::Ph1::SynchronGeneratorVBR::calculateConductanceMatrix() {
	Matrix resistanceMatrix = Matrix::Zero(2,2);
	resistanceMatrix(0,0) = mR_const_1ph.real() + mKa_1ph.real() + mKb_1ph.real();
	resistanceMatrix(0,1) = -mR_const_1ph.imag() - mKa_1ph.imag() + mKb_1ph.imag();
	resistanceMatrix(1,0) = mR_const_1ph.imag() + mKa_1ph.imag() + mKb_1ph.imag();
	resistanceMatrix(1,1) = mR_const_1ph.real() + mKa_1ph.real() - mKb_1ph.real();
	resistanceMatrix = resistanceMatrix * mBase_Z;
	mConductanceMatrix = resistanceMatrix.inverse();
}

void DP::Ph1::SynchronGeneratorVBR::calculateAuxiliarVariables() {	
	mKa = Matrix::Zero(1,3);	
	mKa = mKc * Complex(cos(2. * mThetaMech), sin(2. * mThetaMech));
	mKa_1ph = (mKa * mShiftVector)(0,0);

	mKb = Matrix::Zero(1,3);	
	Real arg = 2. * mThetaMech - 2. * mBase_OmMech * mSimTime ;
	mKb = mKc * Complex(cos(arg), sin(arg));
	mKb_1ph = (mKb * mShiftVectorConj)(0,0);

	mKvbr = Matrix::Zero(1,2);
	mKvbr(0,0) = Complex(cos(mThetaMech - mBase_OmMech * mSimTime), sin(mThetaMech - mBase_OmMech * mSimTime));
	mKvbr(0,1) = -Complex(cos(mThetaMech - mBase_OmMech * mSimTime - PI/2.), sin(mThetaMech - mBase_OmMech * mSimTime - PI/2.));
}

static void setVariableEntries(std::vector<std::pair<UInt, UInt>>& list, Int n, UInt rowIndex, UInt colIndex, Int maxFreq = 1, Int freqIdx = 0) {
	// Assume square matrix
	UInt harmonicOffset = n / maxFreq;
	UInt complexOffset = harmonicOffset / 2;
	UInt harmRow = rowIndex + harmonicOffset * freqIdx;
	UInt harmCol = colIndex + harmonicOffset * freqIdx;

	list.push_back(std::make_pair<UInt, UInt>(rowIndex+0, colIndex+0));
	list.push_back(std::make_pair<UInt, UInt>(rowIndex + complexOffset, colIndex + complexOffset));
	list.push_back(std::make_pair<UInt, UInt>(rowIndex+0, colIndex + complexOffset));
	list.push_back(std::make_pair<UInt, UInt>(rowIndex + complexOffset, colIndex+0));
}

void DP::Ph1::SynchronGeneratorVBR::mnaInitialize(Real omega, 
		Real timeStep, Attribute<Matrix>::Ptr leftVector) {

	Base::ReducedOrderSynchronGenerator<Complex>::mnaInitialize(omega, timeStep, leftVector);

	// TODO: FIX
	Int n = 54;
	
	// upper left
	setVariableEntries(mVariableSystemMatrixEntries, n, mVirtualNodes[0]->matrixNodeIndex(), mVirtualNodes[0]->matrixNodeIndex());

	// bottom right
	setVariableEntries(mVariableSystemMatrixEntries, n, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0));

	// off-diagonal
	setVariableEntries(mVariableSystemMatrixEntries, n, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0, 0));
	setVariableEntries(mVariableSystemMatrixEntries, n, matrixNodeIndex(0, 0), mVirtualNodes[0]->matrixNodeIndex());
	
	mSLog->info("List of index pairs of varying matrix entries: ");
	for (auto indexPair : mVariableSystemMatrixEntries)
		mSLog->info("({}, {})", indexPair.first, indexPair.second);

}

void DP::Ph1::SynchronGeneratorVBR::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	// Stamp voltage source
	Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), mVirtualNodes[1]->matrixNodeIndex(), Complex(-1, 0));
	Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), mVirtualNodes[0]->matrixNodeIndex(), Complex(1, 0));

	// Stamp conductance

	// set upper left block
	Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), mVirtualNodes[0]->matrixNodeIndex(), mConductanceMatrix);

	// set bottom right block
	Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), mConductanceMatrix);

	// Set off-diagonal blocks
	Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0, 0), -mConductanceMatrix);
	Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), mVirtualNodes[0]->matrixNodeIndex(), -mConductanceMatrix);
}

void DP::Ph1::SynchronGeneratorVBR::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	Math::setVectorElement(rightVector, mVirtualNodes[1]->matrixNodeIndex(), mEvbr);
}

void DP::Ph1::SynchronGeneratorVBR::mnaPostStep(const Matrix& leftVector) {
	// update armature voltage and current
	mIntfVoltage(0, 0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	mIntfCurrent(0, 0) = Math::complexFromVectorElement(leftVector, mVirtualNodes[1]->matrixNodeIndex());

	// convert armature voltage into dq reference frame
	Matrix parkTransform = get_parkTransformMatrix();
	MatrixComp Vabc_ = mIntfVoltage(0, 0) * mShiftVector * Complex(cos(mNomOmega * mSimTime), sin(mNomOmega * mSimTime));
	Matrix Vabc = Matrix(3,1);
	Vabc << Vabc_(0,0).real(), Vabc_(1,0).real(), Vabc_(2,0).real();
	mVdq = parkTransform * Vabc / mBase_V_RMS;

	// convert armature current into dq reference frame
	MatrixComp Iabc_ = mIntfCurrent(0, 0) * mShiftVector * Complex(cos(mNomOmega * mSimTime), sin(mNomOmega * mSimTime));
	Matrix Iabc = Matrix(3,1);
	Iabc << Iabc_(0,0).real(), Iabc_(1,0).real(), Iabc_(2,0).real();
	mIdq = parkTransform * Iabc / mBase_I_RMS;
}

Matrix DP::Ph1::SynchronGeneratorVBR::get_parkTransformMatrix() {
	Matrix abcToDq0(2, 3);

	abcToDq0 <<
		2./3.*cos(mThetaMech),  2./3.*cos(mThetaMech - 2.*PI/3.),  2./3.*cos(mThetaMech + 2.*PI/3.),
		-2./3.*sin(mThetaMech), -2./3.*sin(mThetaMech - 2.*PI/3.), -2./3.*sin(mThetaMech + 2.*PI/3.);

	return abcToDq0;
}
