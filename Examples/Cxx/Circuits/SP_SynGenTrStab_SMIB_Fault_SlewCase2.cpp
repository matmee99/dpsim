/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>
#include "../Examples.h"

using namespace DPsim;
using namespace CPS;
using namespace CIM::Examples::Grids::SMIB::SlewCase2;

ScenarioConfig smib;

//Switch to trigger fault at generator terminal
Real SwitchOpen = 1e6;
Real SwitchClosed = 1e6;

void SP_1ph_SynGenTrStab_Fault(String simName, Real timeStep, Real finalTime, bool startFaultEvent, bool endFaultEvent, Real startTimeFault, Real endTimeFault, Real cmdInertia, Real cmdDamping) {
	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components
	auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);

	//Synchronous generator ideal model
	auto genPF = SP::Ph1::SynchronGenerator::make("Generator", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	genPF->setParameters(smib.nomPower, smib.nomPhPhVoltRMS, smib.initActivePower, smib.setPointVoltage*smib.t_ratio, PowerflowBusType::PV);
	genPF->setBaseVoltage(smib.Vnom);
	genPF->modifyPowerFlowBusType(PowerflowBusType::PV);

	//Grid bus as Slack
	auto extnetPF = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::debug);
	extnetPF->setParameters(smib.Vnom);
	extnetPF->setBaseVoltage(smib.Vnom);
	extnetPF->modifyPowerFlowBusType(PowerflowBusType::VD);

	//Line
	auto linePF = SP::Ph1::PiLine::make("PiLine", Logger::Level::debug);
	linePF->setParameters(smib.lineResistance, smib.lineInductance, smib.lineCapacitance, smib.lineConductance);
	linePF->setBaseVoltage(smib.Vnom);

	//Switch
	auto faultPF = CPS::SP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	faultPF->setParameters(SwitchOpen, SwitchClosed);
	faultPF->open();

	// Topology
	genPF->connect({ n1PF });
	faultPF->connect({SP::SimNode::GND, n1PF});
	linePF->connect({ n1PF, n2PF });
	extnetPF->connect({ n2PF });
	auto systemPF = SystemTopology(smib.nomFreq,
			SystemNodeList{n1PF, n2PF},
			SystemComponentList{genPF, linePF,extnetPF, faultPF});

	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->addAttribute("v1", n1PF->attribute("v"));
	loggerPF->addAttribute("v2", n2PF->attribute("v"));

	// Simulation
	Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.doInitFromNodesAndTerminals(false);
	simPF.addLogger(loggerPF);
	simPF.run();

	// ----- Dynamic simulation ------
	String simNameSP = simName + "_SP";
	Logger::setLogDir("logs/"+simNameSP);

	// Nodes
	auto n1SP = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2SP = SimNode<Complex>::make("n2", PhaseType::Single);

	// Components
	auto genSP = SP::Ph1::SynchronGeneratorTrStab::make("SynGen", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	genSP->setStandardParametersPU(smib.nomPower, smib.nomPhPhVoltRMS, smib.nomFreq, smib.Xpd*std::pow(smib.t_ratio,2), cmdInertia*smib.H, smib.Rs, cmdDamping*smib.D );
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower= genPF->getApparentPower();
	genSP->setInitialValues(initApparentPower, smib.initMechPower);
	// genSP->setModelFlags(false, false);

	//Grid bus as Slack
	auto extnetSP = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::debug);
	extnetSP->setParameters(smib.Vnom);
	// Line
	auto lineSP = SP::Ph1::PiLine::make("PiLine", Logger::Level::debug);
	lineSP->setParameters(smib.lineResistance, smib.lineInductance, smib.lineCapacitance, smib.lineConductance);

	//Switch
	auto faultSP = SP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	faultSP->setParameters(SwitchOpen, SwitchClosed);
	faultSP->open();

	// // Variable resistance Switch
	// auto faultSP = SP::Ph1::varResSwitch::make("Br_fault", Logger::Level::debug);
	// faultSP->setParameters(SwitchOpen, SwitchClosed);
	// faultSP->setInitParameters(timeStep);
	// faultSP->open();

	// Topology
	genSP->connect({ n1SP });
	faultSP->connect({SP::SimNode::GND, n1SP });
	lineSP->connect({ n1SP, n2SP });
	extnetSP->connect({ n2SP });
	auto systemSP = SystemTopology(smib.nomFreq,
			SystemNodeList{n1SP, n2SP},
			SystemComponentList{genSP, lineSP, extnetSP, faultSP});

	// Initialization of dynamic topology
	CIM::Reader reader(simNameSP, Logger::Level::debug);
	reader.initDynamicSystemTopologyWithPowerflow(systemPF, systemSP);


	// Logging
	auto loggerSP = DataLogger::make(simNameSP);
	loggerSP->addAttribute("v1", n1SP->attribute("v"));
	loggerSP->addAttribute("v2", n2SP->attribute("v"));
	//gen
	loggerSP->addAttribute("Ep", genSP->attribute("Ep"));
	loggerSP->addAttribute("v_gen", genSP->attribute("v_intf"));
	loggerSP->addAttribute("i_gen", genSP->attribute("i_intf"));
	loggerSP->addAttribute("wr_gen", genSP->attribute("w_r"));
	loggerSP->addAttribute("delta_r_gen", genSP->attribute("delta_r"));
	loggerSP->addAttribute("P_elec", genSP->attribute("P_elec"));
	loggerSP->addAttribute("P_mech", genSP->attribute("P_mech"));
	//Switch
	loggerSP->addAttribute("i_fault", faultSP->attribute("i_intf"));
	//line
	loggerSP->addAttribute("v_line", lineSP->attribute("v_intf"));
	loggerSP->addAttribute("i_line", lineSP->attribute("i_intf"));
	//slack
	// loggerSP->addAttribute("v_slack", extnetSP->attribute("v_intf"));
	// loggerSP->addAttribute("i_slack", extnetSP->attribute("i_intf"));


	Simulation simSP(simNameSP, Logger::Level::debug);
	simSP.setSystem(systemSP);
	simSP.setTimeStep(timeStep);
	simSP.setFinalTime(finalTime);
	simSP.setDomain(Domain::SP);
	simSP.addLogger(loggerSP);
	// simSP.doSystemMatrixRecomputation(true);

	// Events
	if (startFaultEvent){
		auto sw1 = SwitchEvent::make(startTimeFault, faultSP, true);

		simSP.addEvent(sw1);
	}

	if(endFaultEvent){
		std::cout<< endFaultEvent << std::endl;
		// auto sw2 = SwitchEvent::make(endTimeFault, faultSP, false);
		// simSP.addEvent(sw2);

	}

	simSP.run();
}

int main(int argc, char* argv[]) {


	//Simultion parameters
	String simName="SP_SynGenTrStab_SMIB_Fault";
	Real finalTime = 30;
	Real timeStep = 0.001;
	Bool startFaultEvent=true;
	Bool endFaultEvent=true;
	Real startTimeFault=10;
	Real endTimeFault=10.2;
	Real cmdInertia= 1.0;
	Real cmdDamping=1.0;

	CommandLineArgs args(argc, argv);
	if (argc > 1) {
		timeStep = args.timeStep;
		finalTime = args.duration;
		if (args.name != "SPsim")
			simName = args.name;
		if (args.options.find("SCALEINERTIA") != args.options.end())
			cmdInertia = args.options["SCALEINERTIA"];
		if (args.options.find("SCALEDAMPING") != args.options.end())
			cmdDamping = args.options["SCALEDAMPING"];
		if (args.options.find("STARTTIMEFAULT") != args.options.end())
			startTimeFault = args.options["STARTTIMEFAULT"];
		if (args.options.find("ENDTIMEFAULT") != args.options.end())
			endTimeFault = args.options["ENDTIMEFAULT"];
	}

	SP_1ph_SynGenTrStab_Fault(simName, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, cmdInertia, cmdDamping);
}