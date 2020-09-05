/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search.problems.Scheduling;

import org.moeaframework.core.Solution;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.problem.SystemArchitectureProblem;
import seakers.architecture.util.IntegerVariable;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.BaseParams;
import seakers.vassar.problems.Scheduling.Architecture;
import seakers.vassar.problems.Scheduling.SchedulingParams;
import seakers.engineerserver.server.VASSARInterfaceHandler;
import seakers.engineerserver.thriftinterface.*;
import java.util.*;

import java.util.concurrent.ExecutionException;


/**
 * An assigning problem to optimize the allocation of n instruments to m orbits.
 * Also can choose the number of satellites per orbital plane. Objectives are
 * cost and scientific benefit
 *
 * @author nozomihitomi
 */
public class SchedulingProblem extends AbstractProblem implements SystemArchitectureProblem {

    private final String problem;

    private final VASSARInterfaceHandler interfaceHandler;

    private final int numMissions;

    private List<BinaryInputArchitecture> inputArchitectures;

    private List<MissionMeasurements> historicalInfo;

    public SchedulingProblem(String problem, VASSARInterfaceHandler interfaceHandler, int numMissions, List<BinaryInputArchitecture> inputArchitectures, List<MissionMeasurements> historicalInfo) {
        super(numMissions, 2);
        this.problem = problem;
        this.interfaceHandler = interfaceHandler;
        this.numMissions = numMissions;
        this.inputArchitectures = inputArchitectures;
        this.historicalInfo = historicalInfo;
    }

    @Override
    public void evaluate(Solution sltn) {
        SchedulingArchitecture arch = (SchedulingArchitecture) sltn;
        evaluateArch(arch);
        System.out.println(String.format("Arch %s Data Continuity = %10f; Fairness = %10f",
                arch.toString(), arch.getObjective(0), arch.getObjective(1)));
    }

    private void evaluateArch(SchedulingArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {

            int missions = numMissions;
            int[] missionScheduling = new int[missions];

            for (int i = 0; i < missions; i++) {
                missionScheduling[i] = ((IntegerVariable) arch.getVariable(i)).getValue();
            }
            int startYear = 2006;
            List<MissionMeasurements> allMeas = new ArrayList<>();
            for (int i = 0; i < inputArchitectures.size(); i++) {
                List<String> measNames = this.interfaceHandler.getMeasurements(problem, inputArchitectures.get(i));
                List<Double> panelScores = this.interfaceHandler.getPanelScoresForArch(problem, inputArchitectures.get(i));
                MissionMeasurements missionMeas = new MissionMeasurements(measNames,startYear+missionScheduling[i], startYear+missionScheduling[i]+2, panelScores);
                allMeas.add(missionMeas);
            }
            double dataContinuityScore = this.interfaceHandler.evaluateDataContinuityScore(allMeas,historicalInfo);
            double fairnessScore = this.interfaceHandler.evaluateFairnessScore(allMeas);


            arch.setObjective(0, -dataContinuityScore); //negative because MOEAFramework assumes minimization problems
            arch.setObjective(1, -fairnessScore);
            arch.setAlreadyEvaluated(true);
        }
    }

    @Override
    public Solution newSolution() {
        return new SchedulingArchitecture(numMissions, 2);
    }

}
