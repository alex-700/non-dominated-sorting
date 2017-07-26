package ru.ifmo.jmh;

import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import ru.ifmo.JensenFortinBuzdalov;
import ru.ifmo.NonDominatedSortingFactory;

public class JMH_JensenFortinBuzdalov_Fenwick extends AbstractBenchmark {
    @Override
    protected NonDominatedSortingFactory getFactory() {
        return JensenFortinBuzdalov.getFenwickSweepImplementation();
    }

    public static void main(String[] args) throws RunnerException {
        Options opt = new OptionsBuilder()
                .include(JMH_JensenFortinBuzdalov_Fenwick.class.getSimpleName()).build();
        new Runner(opt).run();
    }
}