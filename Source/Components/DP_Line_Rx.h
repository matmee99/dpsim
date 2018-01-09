/** RX Line
 *
 * @file
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
 *
 * DPsim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/

#pragma once

#include "Base.h"

namespace DPsim {
namespace Components {
namespace DP {

	class RxLine : public Components::Base {

	protected:
		Real mResistance;
		Real mConductance;
		Real mVoltageAtNode1Re;
		Real mVoltageAtNode1Im;
		Real mVoltageAtNode2Re;
		Real mVoltageAtNode2Im;

		Real mInductance;
		Real mDeltaVre;
		Real mDeltaVim;
		Real mCurrRe;
		Real mCurrIm;
		Real mCurEqRe;
		Real mCurEqIm;
		Real mGlr;
		Real mGli;
		Real mPrevCurFacRe;
		Real mPrevCurFacIm;

		enum LineTypes {
			RxLine2Node,
			RxLine3Node
		} mType;

		Real correctr, correcti;
		Real cureqr_ind, cureqi_ind;
		Real deltavr_ind;
		Real deltavi_ind;
		Real glr_ind, gli_ind;
		Real currr_ind;
		Real curri_ind;

	public:
		RxLine(String name, Int node1, Int node2, Real resistance, Real inductance, LineTypes type = LineTypes::RxLine3Node);

		void init(SystemModel& system);
		void applySystemMatrixStamp(SystemModel& system);
		void applyRightSideVectorStamp(SystemModel& system) { }
		void step(SystemModel& system, Real time);
		void postStep(SystemModel& system);
		Complex getCurrent(SystemModel& system);
	};
}
}
}