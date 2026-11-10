package scratch.anne.sra.engine.logicTree;

import java.util.Arrays;
import java.util.Map;
import java.util.Objects;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.domain.risk.LossExceedanceCurve;
import scratch.kevin.ucerf3.eal.LossCOV_Model;

public final class LecAccumulation {

    private LecAccumulation() {}

    public static LossExceedanceCurve discreteRank(Map<RuptureKey, Double> branchRatesByRupture,
                                                     Map<RuptureKey, Double> fixedLossLookup) {

        double[] losses = branchRatesByRupture.keySet().stream()
                .map(fixedLossLookup::get)
                .filter(Objects::nonNull)
                .mapToDouble(Double::doubleValue)
                .distinct()
                .sorted()
                .toArray();

        double[] rateByThreshold = new double[losses.length];
        RuptureKey[] keyByThreshold = new RuptureKey[losses.length];

        for (Map.Entry<RuptureKey, Double> entry : branchRatesByRupture.entrySet()) {
            RuptureKey key = entry.getKey();
            double rate = entry.getValue();
            if (rate == 0d)
                continue;

            Double loss = fixedLossLookup.get(key);
            if (loss == null)
                continue;

            int idx = Arrays.binarySearch(losses, loss);
            rateByThreshold[idx] += rate;
            keyByThreshold[idx] = key;
        }

        double[] cumulativeRate = new double[losses.length];
        double running = 0d;
        for (int i = losses.length - 1; i >= 0; i--) {
            running += rateByThreshold[i];
            cumulativeRate[i] = running;
        }

        return new LossExceedanceCurve(losses, cumulativeRate, keyByThreshold);
    }

    public static LossExceedanceCurve logNormalCOV(double[] thresholds,
                                                     Map<RuptureKey, Double> branchRatesByRupture,
                                                     Map<RuptureKey, Double> fixedLossLookup,
                                                     LossCOV_Model covModel) {

        double[] cumulative = new double[thresholds.length];

        DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
        for (double t : thresholds)
            xVals.set(t, 0d);

        for (Map.Entry<RuptureKey, Double> entry : branchRatesByRupture.entrySet()) {
            double rate = entry.getValue();
            if (rate == 0d)
                continue;

            Double loss = fixedLossLookup.get(entry.getKey());
            if (loss == null || loss == 0d)
                continue;

            DiscretizedFunc probs = covModel.calcLossExceedanceProbs(xVals, loss);
            for (int k = 0; k < thresholds.length; k++)
                cumulative[k] += rate * probs.getY(k);
        }

        return new LossExceedanceCurve(thresholds, cumulative, null);
    }
}