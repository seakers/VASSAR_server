/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search.problems.Assigning;

import org.moeaframework.core.Solution;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.problem.SystemArchitectureProblem;
import seakers.vassar.Result;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.local.BaseParams;
import seakers.vassar.problems.Assigning.Architecture;
import seakers.vassar.problems.Assigning.AssigningParams;

import java.util.concurrent.ExecutionException;


/**
 * An assigning problem to optimize the allocation of n instruments to m orbits.
 * Also can choose the number of satellites per orbital plane. Objectives are
 * cost and scientific benefit
 *
 * @author nozomihitomi
 */
public class AssigningProblem extends AbstractProblem implements SystemArchitectureProblem {

    private final int[] alternativesForNumberOfSatellites;

    private final String problem;

    private final ArchitectureEvaluationManager evaluationManager;

    private final BaseParams params;

    private final double dcThreshold = 0.5;

    private final double massThreshold = 3000.0; //[kg]

    private final double packingEffThreshold = 0.4; //[kg]

    /**
     * @param alternativesForNumberOfSatellites
     */
    public AssigningProblem(int[] alternativesForNumberOfSatellites, String problem, ArchitectureEvaluationManager evaluationManager, BaseParams params) {
        //2 decisions for Choosing and Assigning Patterns
        super(1 + params.getNumInstr()*params.getNumOrbits(), 2);
        this.problem = problem;
        this.evaluationManager = evaluationManager;
        this.alternativesForNumberOfSatellites = alternativesForNumberOfSatellites;
        this.params = params;
    }

    @Override
    public void evaluate(Solution sltn) {
        AssigningArchitecture arch = (AssigningArchitecture) sltn;
        evaluateArch(arch);
        System.out.println(String.format("Arch %s Science = %10f; Cost = %10f",
                arch.toString(), arch.getObjective(0), arch.getObjective(1)));
    }

    private void evaluateArch(AssigningArchitecture arch) {
        if (!arch.getAlreadyEvaluated()) {
            StringBuilder bitStringBuilder = new StringBuilder(this.getNumberOfVariables());
            for (int i = 1; i < this.getNumberOfVariables(); ++i) {
                bitStringBuilder.append(arch.getVariable(i).toString());
            }

            AbstractArchitecture arch_old;
            if (problem.equalsIgnoreCase("SMAP") || problem.equalsIgnoreCase("ClimateCentric")) {
                // Generate a new architecture
                arch_old = new Architecture(bitStringBuilder.toString(), 1, (AssigningParams)params);
            }
            else {
                throw new IllegalArgumentException("Unrecorgnizable problem type: " + problem);
            }

            try {
                Result result = this.evaluationManager.evaluateArchitectureAsync(arch_old, "Slow").get();
                arch.setObjective(0, -result.getScience()); //negative because MOEAFramework assumes minimization problems

                double cost = result.getCost();
                arch.setObjective(1, cost); //normalize cost to maximum value
                arch.setAlreadyEvaluated(true);
            }
            catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public Solution newSolution() {
        return new AssigningArchitecture(alternativesForNumberOfSatellites, params.getNumInstr(), params.getNumOrbits(), 2);
    }

}
