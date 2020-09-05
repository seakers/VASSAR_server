package seakers.engineerserver.search.problems.Scheduling;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import seakers.architecture.util.IntegerVariable;
import seakers.vassar.architecture.AbstractArchitecture;
import seakers.vassar.BaseParams;
import seakers.vassar.problems.Scheduling.Architecture;
import seakers.vassar.problems.Scheduling.ArchitectureGenerator;
import seakers.vassar.problems.Scheduling.SchedulingParams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SchedulingInitialization implements Initialization{

    protected final Problem problem;
    protected final int populationSize;
    protected List<Solution> injectedSolutions;
    protected SchedulingParams params;

    private static final boolean debug = false;

    public SchedulingInitialization(Problem problem, int populationSize, BaseParams params){
        this.problem = problem;
        this.populationSize = populationSize;
        this.injectedSolutions = new ArrayList<>();
    }

    public SchedulingInitialization(Problem problem, int populationSize, List<Solution> injectedSolution, BaseParams params){
        this.problem = problem;
        this.populationSize = populationSize;
        this.injectedSolutions = injectedSolution;
    }

    public Solution[] initialize(){

        Solution[] population;

        if (this.populationSize <= this.injectedSolutions.size()) {
            population = this.injectedSolutions.toArray(new Solution[0]);

        }
        else {
            int numRandomlyGeneratedSol = this.populationSize - this.injectedSolutions.size();

            ArchitectureGenerator generator = new ArchitectureGenerator(this.populationSize);
            List<AbstractArchitecture> randomArchs = generator.generateRandomPopulation(numRandomlyGeneratedSol);

            List<Solution> out = new ArrayList<>(injectedSolutions);

            for (AbstractArchitecture arch:randomArchs) {
                Architecture paArch = (Architecture)arch;

                int[] scheduling = paArch.getSchedule();

                SchedulingArchitecture sol = (SchedulingArchitecture) problem.newSolution();
                for (int i = 0; i < scheduling.length; i++) {
                    ((IntegerVariable) sol.getVariable(i)).setValue(scheduling[i]);
                }
                out.add(sol);
            }
            population = out.toArray(new Solution[0]);
        }

        if (debug) {
            for (Solution sol:population){
                int[] scheduling = new int[params.getNumMiss()];
                for (int i = 0; i < scheduling.length; ++i) {
                    scheduling[i] = ((IntegerVariable) sol.getVariable(i)).getValue();
                }
                System.out.println(Arrays.toString(scheduling));
            }
        }

        return  population;
    }
}
