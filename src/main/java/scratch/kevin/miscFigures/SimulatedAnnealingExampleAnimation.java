package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.AnimatedGIFRenderer;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.FullySegmentedRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.InversionState;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.ProgressTrackingCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.Builder;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class SimulatedAnnealingExampleAnimation {

	public static void main(String[] args) throws IOException {
		Location l10 = new Location(34, -118);
		Location l11 = LocationUtils.location(l10, 0d, 120d);
		
		String sectJSON =
				"    {\n"+
				"      \"type\": \"Feature\",\n"+
				"      \"id\": 0,\n"+
				"      \"properties\": {\n"+
				"        \"FaultID\": 0,\n"+
				"        \"FaultName\": \"Test Fault 1\",\n"+
				"        \"DipDeg\": 90.0,\n"+
				"        \"Rake\": 0.0,\n"+
				"        \"LowDepth\": 15.0,\n"+
				"        \"UpDepth\": 0.0,\n"+
				"        \"SlipRate\": 10,\n"+
				"        \"SlipRateStdDev\": 1\n"+
				"      },\n"+
				"      \"geometry\": {\n"+
				"        \"type\": \"LineString\",\n"+
				"        \"coordinates\": [\n"+
				"          [\n"+
				"            "+l10.getLongitude()+",\n"+
				"            "+l10.getLatitude()+"\n"+
				"          ],\n"+
				"          [\n"+
				"            "+l11.getLongitude()+",\n"+
				"            "+l11.getLatitude()+"\n"+
				"          ]\n"+
				"        ]\n"+
				"      }\n"+
				"    }";
		
		File outputDir = new File("/home/kevin/Documents/misc_figures/sa_anim");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(Feature.fromJSON(sectJSON));
//		double subSectLen = sect.getOrigDownDipWidth()*0.5;
//		int minSubSectsPerRup = 2;
		double subSectLen = sect.getOrigDownDipWidth()*0.25;
		int minSubSectsPerRup = 4;
		
		List<GeoJSONFaultSection> subSects = sect.getSubSectionsList(subSectLen, 0, 2);
		System.out.println("Built "+subSects.size()+" sub sects");
		
		RupSetScalingRelationship scale = NSHM23_ScalingRelationships.LOGA_C4p2;
		
		FullySegmentedRupSetConfig config = new RuptureSets.FullySegmentedRupSetConfig(subSects, scale);
		config.setMinSectsPerParent(minSubSectsPerRup);
		FaultSystemRupSet rupSet = config.build(1);
		
		System.out.println("Rupture set has "+rupSet.getNumRuptures()+" ruptures");
		System.out.println("Magnitude range: ["+(float)rupSet.getMinMag()+", "+(float)rupSet.getMaxMag()+"]");
		
		double bVal = 1d;
		
		Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal);
		mfdBuilder.magDepDefaultRelStdDev(M->0.1*Math.max(1, Math.pow(10, bVal*0.5*(M-6))));
		SupraSeisBValInversionTargetMFDs targetMFDs = mfdBuilder.build();
		
		UncertainIncrMagFreqDist totUnboundedTarget = (UncertainIncrMagFreqDist) targetMFDs.getTotalOnFaultSupraSeisMFD();
		UncertainBoundedIncrMagFreqDist totTarget = totUnboundedTarget.estimateBounds(UncertaintyBoundType.ONE_SIGMA);
		
		long iters = 500000;
		int numFrames = 100;
		long subIters = iters/numFrames;
		double fps = 8;
		DecimalFormat frameDF = new DecimalFormat("000");
		
		Range energyXRange = new Range(0d, (double)iters);
		Range energyYRange = new Range(0d, 150);
		
		double mfdMax = totTarget.getUpperMaxY();
		double mfdMin = 1d;
		for (int i=0; i<totTarget.size(); i++) {
			double lowerY = totTarget.getLowerY(i);
			double y = totTarget.getY(i);
			if (lowerY > 0)
				mfdMin = Math.min(mfdMin, lowerY);
			else if (y > 0)
				mfdMin = Math.min(mfdMin, y);
		}
		mfdMax = Math.pow(10, Math.ceil(Math.log10(mfdMax)));
		mfdMin = Math.pow(10, Math.floor(Math.log10(mfdMin)));
		Range mfdXRange = new Range(6.5d, 7.5d);
		Range mfdYRange = new Range(mfdMin, mfdMax);
		
		Range slipXRange = new Range(0d, sect.getTraceLength());
		Range slipYRange = new Range(0d, sect.getOrigAveSlipRate()*1.25);
		
		CompletionCriteria completion = new IterationCompletionCriteria(iters);
		CompletionCriteria subCompletion = new IterationCompletionCriteria(subIters);
		
		List<InversionConstraint> constraints = new ArrayList<>();
		constraints.add(new SlipRateInversionConstraint(1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, rupSet));
		constraints.add(new MFDInversionConstraint(rupSet, 1d, false, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, List.of(totTarget)));
		
		InversionConfiguration.Builder invConfigBuilder = InversionConfiguration.builder(constraints, completion);
		invConfigBuilder.threads(1);
		invConfigBuilder.variablePertubationBasis(Inversions.getDefaultVariablePerturbationBasis(rupSet));
		InversionConfiguration invConvig = invConfigBuilder.build();
		
		InversionInputGenerator inputGen = new InversionInputGenerator(rupSet, invConvig);
		inputGen.generateInputs(true);
		SimulatedAnnealing sa = invConvig.buildSA(inputGen);
		sa.setRandom(new Random(1234l*(long)rupSet.getNumRuptures()));
//		((SerialSimulatedAnnealing)sa).set
		
		ArbitrarilyDiscretizedFunc totEnergyIters = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc slipEnergyIters = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc mfdEnergyIters = new ArbitrarilyDiscretizedFunc();
		totEnergyIters.set(0d, energyYRange.getUpperBound()*100);
		slipEnergyIters.set(0d, energyYRange.getUpperBound()*100);
		mfdEnergyIters.set(0d, energyYRange.getUpperBound()*100);
		
		File framesDir = new File(outputDir, "frames");
		Preconditions.checkState(framesDir.exists() || framesDir.mkdir());
		
		AnimatedGIFRenderer gifRenderer = new AnimatedGIFRenderer(new File(outputDir, "animation.gif"), fps, true);
		
		PlotSpec slipPlot = buildSlipPlot(rupSet, null);
		PlotSpec mfdPlot = buildMFDPlot(rupSet, null, totTarget);
		PlotSpec energyPlot = buildEnergyPlot(totEnergyIters, slipEnergyIters, mfdEnergyIters, 0l, energyXRange, energyYRange);
		
		int frameCount = 0;
		BufferedImage frame = plotCombined(slipPlot, slipXRange, slipYRange, mfdPlot, mfdXRange, mfdYRange,
				energyPlot, energyXRange, energyYRange);
		gifRenderer.writeFrame(frame);
		ImageIO.write(frame, "png", new File(framesDir, "frame_"+frameDF.format(frameCount++)+".png"));
		
		System.out.println("Annealing");
		long iter = 0;
		InversionState state = null;
		while (iter<iters) {
			state = sa.iterate(subCompletion);
			iter += subIters;
			System.out.println("Energy after "+iter+" iters: "+state.energy[0]);
			
			totEnergyIters.set((double)iter, state.energy[0]);
			Preconditions.checkState(state.energy.length == 6);
			slipEnergyIters.set((double)iter, state.energy[4]);
			mfdEnergyIters.set((double)iter, state.energy[5]);
			
			FaultSystemSolution sol = new FaultSystemSolution(rupSet, Arrays.copyOf(sa.getBestSolution(), rupSet.getNumRuptures()));
			
			slipPlot = buildSlipPlot(rupSet, sol);
			mfdPlot = buildMFDPlot(rupSet, sol, totTarget);
			energyPlot = buildEnergyPlot(totEnergyIters, slipEnergyIters, mfdEnergyIters, iter, energyXRange, energyYRange);
			
			frame = plotCombined(slipPlot, slipXRange, slipYRange, mfdPlot, mfdXRange, mfdYRange,
					energyPlot, energyXRange, energyYRange);
			gifRenderer.writeFrame(frame);
			ImageIO.write(frame, "png", new File(framesDir, "frame_"+frameDF.format(frameCount++)+".png"));
		}
		System.out.println("Final energy: "+state.energy[0]);
	}
	
	private static Color TARGET_COLOR = Color.GRAY;
	private static Color SOL_COLOR = Color.BLACK;
	private static float LINE_THICKNESS = 4f;
	private static Color UNCERT_COLOR = new Color(Color.CYAN.getRed(), Color.CYAN.getGreen(), Color.CYAN.getBlue(), 100);
	private static boolean PLOT_UNCERT = false;
	
	private static PlotSpec buildSlipPlot(FaultSystemRupSet rupSet,
			FaultSystemSolution sol) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		SectSlipRates sectSlipRates = rupSet.getSectSlipRates();
		SolutionSlipRates solSlipRates = sol == null ? null : SolutionSlipRates.calc(sol,
				rupSet.requireModule(AveSlipModule.class), rupSet.requireModule(SlipAlongRuptureModel.class));
		XY_DataSet targetSlips = new DefaultXY_DataSet();
		XY_DataSet solSlips = new DefaultXY_DataSet();
		double das = 0d;
		for (int s=0; s<rupSet.getNumSections(); s++) {
			double target = sectSlipRates.getSlipRate(s)*1e3;
			double solution = sol == null ? 0 : solSlipRates.get(s)*1e3;
			
			targetSlips.set(das, target);
			solSlips.set(das, solution);
			das += rupSet.getFaultSectionData(s).getTraceLength();
			targetSlips.set(das, target);
			solSlips.set(das, solution);
		}
		
		targetSlips.setName("Target Slip Rate");
		funcs.add(targetSlips);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, TARGET_COLOR));
		
		solSlips.setName("Solution Slip Rate");
		funcs.add(solSlips);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, SOL_COLOR));
		
		PlotSpec spec = new PlotSpec(funcs, chars,
				"Slip Rates", "Distance Along Strike (km)", "Slip Rate (mm/yr)");
		spec.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
		
		return spec;
	}
	
	private static PlotSpec buildMFDPlot(FaultSystemRupSet rupSet,
			FaultSystemSolution sol,
			UncertainBoundedIncrMagFreqDist targetMFD) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		targetMFD.setName("Target MFD");
		
		PlotCurveCharacterstics targetChar = new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, TARGET_COLOR);
		
		funcs.add(targetMFD);
		chars.add(targetChar);
		
		if (PLOT_UNCERT) {
			// copy with uncertainties
			UncertainBoundedIncrMagFreqDist targetCopy = targetMFD.deepClone();
			targetCopy.setName(targetMFD.getBoundName());
			funcs.add(targetCopy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, UNCERT_COLOR));
			
			// another copy on top of uncertainties
			targetCopy = targetMFD.deepClone();
			targetCopy.setName(null);
			funcs.add(targetMFD);
			chars.add(targetChar);
		}
		
		// solution MFD
		IncrementalMagFreqDist solMFD;
		if (sol == null)
			// empty
			solMFD = new IncrementalMagFreqDist(targetMFD.getMinX(), targetMFD.size(), targetMFD.getDelta());
		else
			solMFD = sol.calcTotalNucleationMFD(targetMFD.getMinX(), targetMFD.getMaxX(), targetMFD.getDelta());
		solMFD.setName("Solution MFD");
		funcs.add(solMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, SOL_COLOR));
		
		PlotSpec spec = new PlotSpec(funcs, chars,
				"Magnitude-Frequency Distribution", "Magnitude", "Annual Frequency");
		spec.setLegendInset(true);
		
		return spec;
	}
	
	private static PlotSpec buildEnergyPlot(
			ArbitrarilyDiscretizedFunc totEnergyIters,
			ArbitrarilyDiscretizedFunc slipEnergyIters,
			ArbitrarilyDiscretizedFunc mfdEnergyIters,
			long iter, Range xRange, Range yRange) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		totEnergyIters.setName("Total");
		funcs.add(totEnergyIters);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, Color.BLACK));
		
		slipEnergyIters.setName("Slip Rates");
		funcs.add(slipEnergyIters);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS/2f, Color.BLUE));
		
		mfdEnergyIters.setName("MFDs");
		funcs.add(mfdEnergyIters);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS/2f, Color.RED));
		
		// clone on top
		totEnergyIters = totEnergyIters.deepClone();
		totEnergyIters.setName(null);
		funcs.add(totEnergyIters);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, LINE_THICKNESS, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars,
				"Inversion Misfits (Energy)", "Iteration Count", "Energy (L^2 Norm of Misfits)");
		spec.setLegendInset(true);
		
		double annX = xRange.getCentralValue();
		double annY = yRange.getLowerBound() + 0.95*yRange.getLength();
		DecimalFormat iterDF = new DecimalFormat("0");
		iterDF.setGroupingSize(3);
		iterDF.setGroupingUsed(true);
		XYTextAnnotation ann = new XYTextAnnotation("Iteration "+iterDF.format(iter), annX, annY);
		ann.setTextAnchor(TextAnchor.TOP_CENTER);
		ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 30));
		spec.addPlotAnnotation(ann);
		
		return spec;
	}
	
	private static BufferedImage plotCombined(PlotSpec slipPlot, Range slipXRange, Range slipYRange,
			PlotSpec mfdPlot, Range mfdXRange, Range mfdYRange,
			PlotSpec energyPlot, Range energyXRange, Range energyYRange) {
		int width = 1800;
		int height = 1000;
		
		int slipX = 0;
		int slipY = 0;
		int slipHeight = height/2;
		int slipWidth = width/2;
		
		int mfdX = width/2;
		int mfdY = 0;
		int mfdHeight = height;
		int mfdWidth = width/2;
		
		int energyX = 0;
		int energyY = height/2;
		int energyHeight = height/2;
		int energyWidth = width/2;
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(slipPlot, false, false, slipXRange, slipYRange);
		gp.getChartPanel().setSize(slipWidth, slipHeight);
		BufferedImage slipBI = gp.getChartPanel().getChart().createBufferedImage(slipWidth, slipHeight);

		gp.drawGraphPanel(mfdPlot, false, true, mfdXRange, mfdYRange);
		gp.getChartPanel().setSize(mfdWidth, mfdHeight);
		BufferedImage mfdBI = gp.getChartPanel().getChart().createBufferedImage(mfdWidth, mfdHeight);
		
		gp.drawGraphPanel(energyPlot, false, false, energyXRange, energyYRange);
		gp.getChartPanel().setSize(energyWidth, energyHeight);
		BufferedImage energyBI = gp.getChartPanel().getChart().createBufferedImage(energyWidth, energyHeight);
		
		BufferedImage stitched = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		for (int x=0; x<slipWidth; x++)
			for (int y=0; y<slipHeight; y++)
				stitched.setRGB(x+slipX, y+slipY, slipBI.getRGB(x, y));
		for (int x=0; x<mfdWidth; x++)
			for (int y=0; y<mfdHeight; y++)
				stitched.setRGB(x+mfdX, y+mfdY, mfdBI.getRGB(x, y));
		for (int x=0; x<energyWidth; x++)
			for (int y=0; y<energyHeight; y++)
				stitched.setRGB(x+energyX, y+energyY, energyBI.getRGB(x, y));
		
		return stitched;
	}

}
