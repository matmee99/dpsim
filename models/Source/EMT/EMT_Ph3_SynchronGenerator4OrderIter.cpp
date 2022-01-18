/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph3_SynchronGenerator4OrderIter.h>

using namespace CPS;

EMT::Ph3::SynchronGenerator4OrderIter::SynchronGenerator4OrderIter
    (String uid, String name, Logger::Level logLevel)
	: Base::SimpSynchronousGenerator<Real>(uid, name, logLevel) {

	mPhaseType = PhaseType::ABC;
	setTerminalNumber(1);

	// model flags
	mVoltageForm = false;
	
	// model variables
	mIntfVoltage = Matrix::Zero(3,1);
	mIntfCurrent = Matrix::Zero(3,1);
	mEdq0_t = Matrix::Zero(3,1);
	mEdq0_t_pred = Matrix::Zero(3,1);
	mEdq0_t_corr = Matrix::Zero(3,1);
	mdEdq0_t = Matrix::Zero(3,1);
	mdEdq0_t_corr = Matrix::Zero(3,1);

    // Register attributes
	addAttribute<Int>("NIterations", &mNumIterations2 , Flags::read);
	addAttribute<Matrix>("Edq0_t", &mEdq0_t, Flags::read);
}

EMT::Ph3::SynchronGenerator4OrderIter::SynchronGenerator4OrderIter
	(String name, Logger::Level logLevel)
	: SynchronGenerator4OrderIter(name, name, logLevel) {
}

SimPowerComp<Real>::Ptr EMT::Ph3::SynchronGenerator4OrderIter::clone(String name) {
	auto copy = SynchronGenerator4OrderIter::make(name, mLogLevel);
	
	return copy;
}

void EMT::Ph3::SynchronGenerator4OrderIter::specificInitialization() {

	// Initialize matrix of state representation
	mA = Matrix::Zero(3,3);
	mB = Matrix::Zero(3,3);
	mC = Matrix::Zero(3,1);
	calculateStateMatrix();

	// initial voltage behind the transient reactance in the dq0 reference frame
	mEdq0_t(0,0) = mVdq0(0,0) - mIdq0(1,0) * mLq_t;
	mEdq0_t(1,0) = mVdq0(1,0) + mIdq0(0,0) * mLd_t;

	mSLog->info(
		"\n--- Model specific initialization  ---"
		"\nInitial Ed_t (per unit): {:f}"
		"\nInitial Eq_t (per unit): {:f}"
		"\nMax number of iterations: {:d}"
		"\nTolerance: {:f}"
		"\n--- Model specific initialization finished ---",

		mEdq0_t(0,0),
		mEdq0_t(1,0),
		mMaxIter,
		mTolerance
	);
	mSLog->flush();
}

void EMT::Ph3::SynchronGenerator4OrderIter::calculateStateMatrix() {
	if (mVoltageForm) {
		Real Td_t =  mTd0_t * (mLd_t / mLd);
		Real Tq_t =  mTq0_t * (mLq_t / mLq);
		mA << -1. / Tq_t   ,          0.0,		0.0,
	               0.0     ,      -1 / Td_t,	0.0,
				   0.0     ,          0.0,		0.0;
		mB << (1. / Tq_t) * (mLq-mLq_t) / mLq,   				0.0, 					0.0,
		  					0.0, 					(1. / Td_t) * (mLd-mLd_t) / mLd,	0.0,
							0.0, 								0.0,					0.0;
		mC <<               0.0,
		   		(1. / Td_t) * mEf * (mLd_t / mLd),
				   			0.0;
	}
	else {	// Currents form
		mA << -1./mTq0_t	,    	0.0,		0.0,		
		     	   0.0   	,    -1/mTd0_t,		0.0,
				   0.0		,		0.0,		0.0;
		mB <<             0.0                  , (1. / mTq0_t) * (mLq-mLq_t),		0.0,
		  		(-1. / mTd0_t) * (mLd-mLd_t)   ,            0.0             ,		0.0,
				  		  0.0				   ,			0.0				,		0.0;	
		mC <<          0.0,
		 		 (1./mTd0_t) * mEf,
				  	   0.0;
	}
}

void EMT::Ph3::SynchronGenerator4OrderIter::stepInPerUnit() {
	// set number of iteratios equal to zero
	mNumIter = 0;

	// Predictor step (euler)

	//predict mechanical variables
	mElecTorque = mVdq0(0,0) * mIdq0(1,0) + mVdq0(1,0) * mIdq0(0,0);
	mdOmMech = 1 / (2.* mH) * (mMechTorque - mElecTorque);
	mOmMech_pred = mOmMech + mTimeStep * mdOmMech;
	mdDelta = (mOmMech_pred - 1.) * mBase_OmMech;
	mDelta_pred = mDelta + mTimeStep * mdDelta;
	mThetaMech_pred = mThetaMech + mTimeStep * mOmMech_pred * mBase_OmMech;

	//predict voltage behind transient reactance
	if (mVoltageForm)
		mdEdq0_t = mA * mEdq0_t + mB * mVdq0 + mC;
	else
		mdEdq0_t = mA * mEdq0_t + mB * mIdq0 + mC;
	mEdq0_t_pred = mEdq0_t + mTimeStep * mdEdq0_t;

	// predict armature currents for at t=k+1
	mIdq0(0,0) = (mEdq0_t_pred(1,0) - mVdq0(1,0) ) / mLd_t;
	mIdq0(1,0) = (mVdq0(0,0) - mEdq0_t_pred(0,0) ) / mLq_t;
	mIdq0(2,0) = 0.0;

	// convert currents into the abc domain
	mIntfCurrent = inverseParkTransform(mThetaMech_pred, mIdq0);
	mIntfCurrent = mIntfCurrent * mBase_I;
}

void EMT::Ph3::SynchronGenerator4OrderIter::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	Math::setVectorElement(rightVector, matrixNodeIndex(0,0), mIntfCurrent(0, 0));
	Math::setVectorElement(rightVector, matrixNodeIndex(0,1), mIntfCurrent(1, 0));
	Math::setVectorElement(rightVector, matrixNodeIndex(0,2), mIntfCurrent(2, 0));
}

void EMT::Ph3::SynchronGenerator4OrderIter::correctorStep() {
	// corrector step (trapezoidal rule)

	if (mNumIter==0)
		return;

	//predict mechanical variables
	mElecTorque_corr = mVdq0(0,0) * mIdq0(1,0) + mVdq0(1,0) * mIdq0(0,0);
	mdOmMech_corr = 1 / (2.* mH) * (mMechTorque - mElecTorque_corr);
	mOmMech_corr = mOmMech + mTimeStep / 2. * (mdOmMech + mdOmMech_corr);
	mdDelta_corr = (mOmMech_corr - 1.) * mBase_OmMech;
	mDelta_corr = mDelta + mTimeStep / 2. * (mdDelta + mdDelta_corr);
	mThetaMech_corr = mThetaMech + mTimeStep / 2. *(mOmMech + mOmMech_corr) * mBase_OmMech;

	//predict voltage behind transient reactance
	if (mVoltageForm)
		mdEdq0_t_corr = mA * mEdq0_t + mB * mVdq0 + mC;
	else
		mdEdq0_t_corr = mA * mEdq0_t + mB * mIdq0 + mC;
	mEdq0_t_corr = mEdq0_t + mTimeStep / 2 * (mdEdq0_t + mdEdq0_t_corr);

	// armature currents for at t=k+1
	mIdq0(0,0) = (mEdq0_t_corr(1,0) - mVdq0(1,0) ) / mLd_t;
	mIdq0(1,0) = (mVdq0(0,0) - mEdq0_t_corr(0,0) ) / mLq_t;

	// convert currents into the abc domain
	mIntfCurrent = inverseParkTransform(mThetaMech_corr, mIdq0);
	mIntfCurrent = mIntfCurrent * mBase_I;

	// stamp currents
	mnaApplyRightSideVectorStamp(mRightVector);
}

void EMT::Ph3::SynchronGenerator4OrderIter::updateVoltage(const Matrix& leftVector) {	
	mVdq0_prev = mVdq0;
	
	mIntfVoltage(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	mIntfVoltage(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
	mIntfVoltage(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));

	// convert Vdq into abc domain
	mVdq0 = parkTransform(mThetaMech_pred, mIntfVoltage);
	mVdq0 = mVdq0 / mBase_V;
}

bool EMT::Ph3::SynchronGenerator4OrderIter::checkVoltageDifference() {
	if (mNumIter==0) {
		// if no corrector step has been performed yet
		mNumIter = 1;
		return true;
	}
	
	Matrix voltageDifference = mVdq0 - mVdq0_prev;
	if (Math::abs(voltageDifference(0,0)) > mTolerance || Math::abs(voltageDifference(1,0)) > mTolerance) {
		if (mNumIter > mMaxIter) {
			return false;
		} else {
			mNumIter = mNumIter + 1;
			return true;
		}
	} else
		return false;
}

void EMT::Ph3::SynchronGenerator4OrderIter::mnaPostStep(const Matrix& leftVector) {
	// update variables
	mEdq0_t = mEdq0_t_corr;
	mOmMech = mOmMech_corr;
	mThetaMech = mThetaMech_corr;
	mDelta = mDelta_corr;

	//
	mNumIterations2 = mNumIter;
}

Matrix EMT::Ph3::SynchronGenerator4OrderIter::parkTransform(Real theta, const Matrix& abcVector) {
	Matrix dq0Vector(3, 1);
	Matrix abcToDq0(3, 3);

	// Park transform according to Kundur
	abcToDq0 <<
		 2./3.*cos(theta),	2./3.*cos(theta - 2.*PI/3.),  2./3.*cos(theta + 2.*PI/3.),
		-2./3.*sin(theta), -2./3.*sin(theta - 2.*PI/3.), -2./3.*sin(theta + 2.*PI/3.),
		 1./3., 			1./3., 						  1./3.;

	dq0Vector = abcToDq0 * abcVector;

	return dq0Vector;
}

Matrix EMT::Ph3::SynchronGenerator4OrderIter::inverseParkTransform(Real theta, const Matrix& dq0Vector) {
	Matrix abcVector(3, 1);
	Matrix dq0ToAbc(3, 3);

	// Park transform according to Kundur
	dq0ToAbc <<
		cos(theta), 		   -sin(theta), 		   1.,
		cos(theta - 2.*PI/3.), -sin(theta - 2.*PI/3.), 1.,
		cos(theta + 2.*PI/3.), -sin(theta + 2.*PI/3.), 1.;

	abcVector = dq0ToAbc * dq0Vector;

	return abcVector;
}
