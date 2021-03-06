package ru.ifmo.nds.jmh;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;
import ru.ifmo.nds.IdCollection;
import ru.ifmo.nds.NonDominatedSorting;

import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.SECONDS)
@Timeout(time = 1, timeUnit = TimeUnit.HOURS)
@Warmup(time = 6, iterations = 1)
@Measurement(time = 1, iterations = 1)
@Fork(5)
public class UniformCorrelated {
    private static final int INSTANCES = 10;

    @Param("The algorithm should be set explicitly")
    private String algorithmId;
    private NonDominatedSorting sorting;
    private double[][][] dataset;
    private int[] ranks;

    @Param({"10", "100", "1000", "10000"})
    private int n;

    @Param({"2", "3", "5", "10"})
    private int d;

    @Setup
    public void initializeSorterAndData() {
        sorting = IdCollection.getNonDominatedSortingFactory(algorithmId).getInstance(n, d);
        ranks = new int[n];
        dataset = new double[INSTANCES][n][d];
        Random random = new Random(Arrays.hashCode(new int[] {n, d}));
        for (int i = 0; i < INSTANCES; ++i) {
            fill(random, dataset[i], i % 2 == 0 ? 1 : d - 2);
        }
    }

    private void fill(Random random, double[][] instance, int x) {
        for (int i = 0; i < n; ++i) {
            double first = random.nextDouble();
            for (int k = 0; k < d; ++k) {
                instance[i][k] = k == x ? -first : first;
            }
        }
        Collections.shuffle(Arrays.asList(instance), random);
    }

    @OperationsPerInvocation(INSTANCES)
    @Benchmark
    public void benchmarkCall(Blackhole blackhole) {
        for (double[][] instance : dataset) {
            Arrays.fill(ranks, 0);
            sorting.sort(instance, ranks);
            blackhole.consume(ranks);
        }
    }

    @TearDown
    public void destroySorter() {
        sorting.close();
    }
}
