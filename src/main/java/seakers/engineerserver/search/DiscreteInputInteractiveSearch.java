/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seakers.engineerserver.search;


import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import io.lettuce.core.RedisClient;
import io.lettuce.core.api.StatefulRedisConnection;
import io.lettuce.core.api.sync.RedisCommands;
import io.lettuce.core.pubsub.StatefulRedisPubSubConnection;
import io.lettuce.core.pubsub.api.sync.RedisPubSubCommands;
import seakers.engineerserver.javaInterface.DiscreteInputArchitecture;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;
import seakers.architecture.util.IntegerVariable;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 *
 * @author nozomihitomi
 */
public class DiscreteInputInteractiveSearch implements Callable<Algorithm> {

    private final Algorithm alg;
    private final TypedProperties properties;
    private final String username;
    private final RedisClient redisClient;
    private boolean isStopped;
    private ConcurrentLinkedQueue<Integer> messageQueue;

    public DiscreteInputInteractiveSearch(Algorithm alg, TypedProperties properties, String username, RedisClient redisClient, ConcurrentLinkedQueue<Integer> messageQueue) {
        this.alg = alg;
        this.properties = properties;
        this.username = username;
        this.redisClient = redisClient;
        this.isStopped = false;
        this.messageQueue = messageQueue;
    }

    @Override
    public Algorithm call() {

        int populationSize = (int) properties.getDouble("populationSize", 600);
        int maxEvaluations = (int) properties.getDouble("maxEvaluations", 10000);

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations) && !isStopped) {
            Integer message = messageQueue.poll();
            if (message != null) {
                if (message == 0) {
                    this.isStopped = true;
                }
            }
            alg.step();
            Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
            StatefulRedisConnection<String, String> connection = redisClient.connect();
            RedisCommands<String, String> syncCommands = connection.sync();

            for(int i = 0; i < pop.size(); i++){
                Solution s = pop.get(i);
                s.setAttribute("NFE", alg.getNumberOfEvaluations());
                // Send the new architectures through REDIS
                // But first, turn it into something easier in JSON
                DiscreteInputArchitecture json_arch = new DiscreteInputArchitecture();
                json_arch.inputs = new ArrayList<>();
                json_arch.outputs = new ArrayList<>();
                for (int j = 0; j < s.getNumberOfVariables(); ++j) {
                    IntegerVariable var = (IntegerVariable)s.getVariable(j);
                    int val = var.getValue();
                    json_arch.inputs.add(val);
                }
                json_arch.outputs.add(-s.getObjective(0));
                json_arch.outputs.add(s.getObjective(1));
                Gson gson = new GsonBuilder().create();
                Long retval = syncCommands.rpush(this.username, gson.toJson(json_arch));
            }
            syncCommands.ltrim(this.username, -1000, -1);
            connection.close();
            // Notify listeners of new architectures in username channel
            StatefulRedisPubSubConnection<String, String> pubsubConnection = redisClient.connectPubSub();
            RedisPubSubCommands<String, String> sync = pubsubConnection.sync();
            sync.publish(this.username, "new_arch");
            pubsubConnection.close();
        }

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        return alg;
    }
}