package search;


import javaInterface.DiscreteInputArchitecture;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Utils {

    public static boolean dominates(List<Double> objectives1, List<Double> objectives2){

        // Assumes Smaller-is-better
        boolean at_least_as_good_as = true;
        boolean better_than_in_one = false;

        for(int i = 0; i < objectives1.size(); i++){

            if(objectives1.get(i) < objectives2.get(i)){
                // First better than Second
                better_than_in_one=true;

            }else if(objectives1.get(i) <= objectives2.get(i)){
                // First is worse than Second
                at_least_as_good_as = false;
            }
        }

        return at_least_as_good_as && better_than_in_one; // First dominates Second
    }

    public static boolean dominates(double[] objectives1, double[] objectives2){

        // Assumes Smaller-is-better
        boolean at_least_as_good_as = true;
        boolean better_than_in_one = false;

        for(int i = 0; i < objectives1.length; i++){

            if(objectives1[i] < objectives2[i]){
                // First better than Second
                better_than_in_one=true;

            }else if(objectives1[i] <= objectives2[i]){
                // First is worse than Second
                at_least_as_good_as = false;
            }
        }

        return at_least_as_good_as && better_than_in_one; // First dominates Second
    }


    public static List<DiscreteInputArchitecture> getFuzzyParetoFrontDiscreteInput(List<DiscreteInputArchitecture> population,
                                                                                          List<Boolean> SIB,
                                                                                          int paretoRank){

        List<DiscreteInputArchitecture> fuzzy_pareto_front = new ArrayList<>();

        ArrayList<DiscreteInputArchitecture> current_population = new ArrayList<>();
        for (DiscreteInputArchitecture s:population) {
            current_population.add(s);
        }

        int iteration=0;
        while(iteration <= paretoRank){

            ArrayList<DiscreteInputArchitecture> features_next_iter = new ArrayList<>();
            for(int i = 0; i < current_population.size(); i++){

                boolean dominated = false;
                DiscreteInputArchitecture f1 = current_population.get(i);
                List<Double> obj1 = new ArrayList<>();
                for(int k = 0; k < f1.getOutputsSize(); k++){
                    if(SIB.get(k)){
                        obj1.add(f1.getOutputs().get(k));
                    }else{
                        obj1.add(-f1.getOutputs().get(k));
                    }
                }

                for(int j = 0; j < current_population.size(); j++){

                    if(i==j) continue;

                    DiscreteInputArchitecture f2 = current_population.get(j);
                    List<Double> obj2 = new ArrayList<>();
                    for(int k = 0; k < f2.getOutputsSize(); k++){
                        if(SIB.get(k)){
                            obj2.add(f1.getOutputs().get(k));
                        }else{
                            obj2.add(-f1.getOutputs().get(k));
                        }
                    }
                    if(dominates(obj2, obj1)){ // f1 is dominated
                        dominated=true;
                        break;
                    }
                }

                if(!dominated){
                    fuzzy_pareto_front.add(f1);
                }else{
                    features_next_iter.add(f1);
                }

            }

            current_population = features_next_iter;
            iteration++;
        }

        return fuzzy_pareto_front;
    }
}
