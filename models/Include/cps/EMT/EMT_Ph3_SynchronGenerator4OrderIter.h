/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cps/Solver/MNASyncGenInterface.h>
#include <cps/Base/Base_SimpSynchronousGenerator.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
	/// @brief 4 Order Synchronous generator model for transient stability analysis
	///
	/// This model is based on Eremia section 2.1.6.
	class SynchronGenerator4OrderIter :
		public Base::SimpSynchronousGenerator<Real>,
		public MNASyncGenInterface,
		public SharedFactory<SynchronGenerator4OrderIter> {
	protected:
		///
		Int mNumIterations2;

		/// sim flags
		bool mVoltageForm;

		// #### Model specific variables ####
		///
		Matrix mVdq0_prev;
		/// previous voltage behind the transient impedance (p.u.)
		Matrix mEdq0_t;
		Matrix mEdq0_t_pred;
		Matrix mEdq0_t_corr;
		/// derivative voltage behind the transient impedance (p.u.)
		Matrix mdEdq0_t;
		Matrix mdEdq0_t_corr;
		///
		Real mElecTorque_corr;
		///
		Real mdOmMech;
		Real mdOmMech_corr;
		Real mOmMech_pred;
		Real mOmMech_corr;
		/// prediction of mechanical system angle
		Real mThetaMech_pred;
		Real mThetaMech_corr;
		///
		Real mDelta_pred;
		Real mDelta_corr;
		Real mdDelta;
		Real mdDelta_corr;

		/// State Matrix x(k+1) = Ax(k) + Bu(k) + C
		/// A Matrix
		Matrix mA;
		/// B Matrix
		Matrix mB;
		/// Constant Matrix
		Matrix mC;

	public:
		///
		SynchronGenerator4OrderIter(String uid, String name, Logger::Level logLevel = Logger::Level::off);
		///
		SynchronGenerator4OrderIter(String name, Logger::Level logLevel = Logger::Level::off);
		///
		SimPowerComp<Real>::Ptr clone(String name);

		// #### General Functions ####
		///
		void specificInitialization();
		///
		void calculateStateMatrix();
		///
		void stepInPerUnit();
		// 
		void correctorStep();
		/// 
		void updateVoltage(const Matrix& leftVector);
		///
		bool checkVoltageDifference();
		///
		Matrix parkTransform(Real theta, const Matrix& abcVector);
		///
		Matrix inverseParkTransform(Real theta, const Matrix& dq0Vector);


		/// Setters
		///
		void useVoltageForm(bool state) {mVoltageForm = state;}	

		// #### MNA Functions ####		
		void mnaApplyRightSideVectorStamp(Matrix& rightVector);
		void mnaPostStep(const Matrix& leftVector);
		void mnaApplySystemMatrixStamp(Matrix& systemMatrix){};
	};
}
}
}
