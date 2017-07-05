package scratch.kevin.ucerf3;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.InversionConfiguration;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.InversionInputGenerator;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.simulatedAnnealing.ThreadedSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

import com.google.common.base.Stopwatch;

public class GRMemTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			InversionModels im = InversionModels.GR_CONSTRAINED;
//			FaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(FaultModels.FM3_1,
//					DeformationModels.GEOLOGIC_PLUS_ABM, im);
			LaughTestFilter filter = LaughTestFilter.getDefault();
			filter.setMaxCmlAzimuthChange(180);
			InversionFaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(filter,
					InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE, FaultModels.FM3_1,
					DeformationModels.GEOLOGIC_PLUS_ABM, ScalingRelationships.ELLSWORTH_B, SlipAlongRuptureModels.TAPERED,
					im);
			InversionConfiguration config = InversionConfiguration.forModel(im, rupSet);
			
			ArrayList<PaleoRateConstraint> paleoConstraints =
					UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList());
			PaleoProbabilityModel paleoModel = InversionInputGenerator.loadDefaultPaleoProbabilityModel();
			
			InversionInputGenerator gen =
					new InversionInputGenerator(rupSet, config, paleoConstraints, null, paleoModel);
			
			System.out.println("Generating inputs...");
			Stopwatch watch = Stopwatch.createStarted();
			if (args.length > 0)
				gen.generateInputs(SparseRCDoubleMatrix2D.class);
			else
				gen.generateInputs();
			System.out.println("Done generating inputs after "
				+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d / 60d)+" mins.");
			
			System.gc();
			
			watch.reset();
			watch.start();
			System.out.println("Column compressing...");
			gen.columnCompress();
			
			System.out.println("Done column compressing after "
					+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d / 60d)+" mins.");
			
			System.gc();
			
			System.out.println("Building annealer!");
			
			ThreadedSimulatedAnnealing tsa = new ThreadedSimulatedAnnealing(gen.getA(), gen.getD(), gen.getInitial(), 0,
					gen.getA_ineq(), gen.getD_ineq(), gen.getMinimumRuptureRates(), 23, TimeCompletionCriteria.getInSeconds(1));
			
			System.out.println("Starting annealing!");
			tsa.iterate(TimeCompletionCriteria.getInMinutes(5));
		} catch (Exception e) {
			System.out.println("Ah damnit - it didn't work. :-(");
			e.printStackTrace();
			System.exit(1);
		}
		System.out.println("Done! Yay!");
		System.exit(0);
	}

}
