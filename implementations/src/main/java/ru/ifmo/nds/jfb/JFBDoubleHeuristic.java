package ru.ifmo.nds.jfb;

import ru.ifmo.nds.util.ArrayHelper;
import ru.ifmo.nds.util.ArraySorter;
import ru.ifmo.nds.util.RankQueryStructureDouble;

import java.util.Arrays;

public class JFBDoubleHeuristic extends JFBDouble {
    public JFBDoubleHeuristic(RankQueryStructureDouble rankQueryStructure, int maximumDimension, int allowedThreads, HybridAlgorithmWrapper hybridWrapper) {
        super(rankQueryStructure, maximumDimension, allowedThreads, hybridWrapper);
    }

    @Override
    public String getName() {
        return super.getName() + "+ heuristic";
    }

    private int heuristic(int n, int dim) {
        boolean[] hasRank0 = new boolean[n];
        hasRank0[0] = true;
        for (int c = 1; c < dim; c++) {
            double cur_min = transposedPoints[c][0];
            for (int i = 1; i < n; i++) {
                double cur = transposedPoints[c][i];
                if (cur_min > cur) {
                    cur_min = cur;
                    hasRank0[i] = true;
                }
            }
        }
        int tail = 0;
        for (int i = 0; i < n; i++) {
            if (hasRank0[i]) {
                indices[tail++] = i;
            }
        }
        int ans = tail;
        for (int i = 0; i < n; i++) {
            if (!hasRank0[i]) {
                indices[tail++] = i;
            }
        }
        return ans;
    }

    @Override
    protected void sortChecked(double[][] points, int[] ranks, int maximalMeaningfulRank) {
        final int n = points.length;
        final int dim = points[0].length;
        Arrays.fill(ranks, 0);
        ArrayHelper.fillIdentity(indices, n);
        sorter.lexicographicalSort(points, indices, 0, n, dim);

        this.maximalMeaningfulRank = maximalMeaningfulRank;

        if (dim == 2) {
            // 2: Special case: binary search.
            twoDimensionalCase(points, ranks);
        } else {
            // 3: General case.
            // 3.1: Moving points in a sorted order to internal structures
            final int newN = ArraySorter.retainUniquePoints(points, indices, this.points, ranks);
            Arrays.fill(this.ranks, 0, newN, 0);

            // 3.2: Transposing points. This should fit in cache for reasonable dimensions.
            for (int i = 0; i < newN; ++i) {
                for (int j = 0; j < dim; ++j) {
                    transposedPoints[j][i] = this.points[i][j];
                }
            }

            postTransposePointHook(newN);
            ArrayHelper.fillIdentity(indices, newN);

            // 3.3: Calling the actual sorting
            int split = heuristic(newN, dim);
            int newN2 = helperB(0, split, split, newN, dim - 1, 0);
            helperA(split, newN2, dim - 1);

            // 3.4: Applying the results back. After that, the argument "ranks" array stops being abused.
            for (int i = 0; i < n; ++i) {
                ranks[i] = this.ranks[ranks[i]];
                this.points[i] = null;
            }
        }
    }
}
