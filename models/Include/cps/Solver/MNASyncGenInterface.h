/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cps/AttributeList.h>
#include <cps/Config.h>
#include <cps/Definitions.h>

namespace CPS {
	/// Interface to be used by synchronous generators
	class MNASyncGenInterface {
	public:
		typedef std::shared_ptr<MNASyncGenInterface> Ptr;
		typedef std::vector<Ptr> List;

		/// Returns true if it needs to iterate
		virtual void correctorStep()=0;
		///
		virtual void updateVoltage(const Matrix& leftVector)=0;
		///
		virtual bool checkVoltageDifference() {return false;}

		/// Setters
		/// Default = 100000
		void setMaxIterations(Int maxIterations) {mMaxIter = maxIterations;}
		/// Default=1e-6
		void setTolerance(Real Tolerance) {mTolerance = Tolerance;}	

		///
		Int mNumIter = 0;
		///
		Int mMaxIter = 10;
		///
		Real mTolerance = 1e-6;
	};
}
