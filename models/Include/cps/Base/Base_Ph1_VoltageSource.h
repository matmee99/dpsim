/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cps/Definitions.h>
#include <cps/AttributeList.h>

namespace CPS {
namespace Base {
namespace Ph1 {
	class VoltageSource {
	public:
		/// Voltage set point [V]
		Attribute<Complex>::Ptr mVoltageRef;
		/// Source frequency [Hz]
		Attribute<Real>::Ptr mSrcFreq;
		/// Sets model specific parameters
		void setParameters(Complex voltageRef, Real srcFreq = -1) {
			**mVoltageRef = voltageRef;
			**mSrcFreq = srcFreq;
		}
	};
}
}
}
