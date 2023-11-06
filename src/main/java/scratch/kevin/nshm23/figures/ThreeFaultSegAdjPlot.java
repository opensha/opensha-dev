package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.poi.hssf.util.HSSFColor.BLACK;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectAreas;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;

public class ThreeFaultSegAdjPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/seg_adjustments");
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		// create fake rupture set

		double supraB = 0.5;
		Shaw07JumpDistProb segModel = Shaw07JumpDistProb.forHorzOffset(1d, 3, 2);

		List<FaultSection> sects = new ArrayList<>();
		double len1 = 50d;
		double len2 = 50d;
		double len3 = 50d;

		double az1 = 0d;
		double az2 = 0d;
		double az3 = 0d;

		//	double dist2 = 10d;
		//	double dist3 = 5d;
		double dist2 = segModel.calcJumpDistance(0.2);
		double dist3 = segModel.calcJumpDistance(0.1);

		Location l10 = new Location(34, -118);
		Location l11 = LocationUtils.location(l10, az1, len1);

		Location l20 = LocationUtils.location(l11, az2, dist2);
		Location l21 = LocationUtils.location(l20, az2, len2);

		Location l30 = LocationUtils.location(l21, az3, dist3);
		Location l31 = LocationUtils.location(l30, az3, len3);

		//	Location l3 = new Location
		String sect1JSON =
				"    {\n"+
						"      \"type\": \"Feature\",\n"+
						"      \"id\": 0,\n"+
						"      \"properties\": {\n"+
						"        \"FaultID\": 0,\n"+
						"        \"FaultName\": \"Test Fault 1\",\n"+
						"        \"DipDeg\": 90.0,\n"+
						"        \"Rake\": 0.0,\n"+
						"        \"LowDepth\": 10.0,\n"+
						"        \"UpDepth\": 0.0,\n"+
						"        \"SlipRate\": 10\n"+
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
		//	System.out.println(sect1JSON);
		sects.add(GeoJSONFaultSection.fromFeature(Feature.fromJSON(sect1JSON)));
		String sect2JSON =
				"    {\n"+
						"      \"type\": \"Feature\",\n"+
						"      \"id\": 1,\n"+
						"      \"properties\": {\n"+
						"        \"FaultID\": 0,\n"+
						"        \"FaultName\": \"Test Fault 2\",\n"+
						"        \"DipDeg\": 90.0,\n"+
						"        \"Rake\": 0.0,\n"+
						"        \"LowDepth\": 10.0,\n"+
						"        \"UpDepth\": 0.0,\n"+
						"        \"SlipRate\": 10\n"+
						"      },\n"+
						"      \"geometry\": {\n"+
						"        \"type\": \"LineString\",\n"+
						"        \"coordinates\": [\n"+
						"          [\n"+
						"            "+l20.getLongitude()+",\n"+
						"            "+l20.getLatitude()+"\n"+
						"          ],\n"+
						"          [\n"+
						"            "+l21.getLongitude()+",\n"+
						"            "+l21.getLatitude()+"\n"+
						"          ]\n"+
						"        ]\n"+
						"      }\n"+
						"    }";
		sects.add(GeoJSONFaultSection.fromFeature(Feature.fromJSON(sect2JSON)));
		String sect3JSON =
				"    {\n"+
						"      \"type\": \"Feature\",\n"+
						"      \"id\": 2,\n"+
						"      \"properties\": {\n"+
						"        \"FaultID\": 0,\n"+
						"        \"FaultName\": \"Test Fault 3\",\n"+
						"        \"DipDeg\": 90.0,\n"+
						"        \"Rake\": 0.0,\n"+
						"        \"LowDepth\": 10.0,\n"+
						"        \"UpDepth\": 0.0,\n"+
						"        \"SlipRate\": 10\n"+
						"      },\n"+
						"      \"geometry\": {\n"+
						"        \"type\": \"LineString\",\n"+
						"        \"coordinates\": [\n"+
						"          [\n"+
						"            "+l30.getLongitude()+",\n"+
						"            "+l30.getLatitude()+"\n"+
						"          ],\n"+
						"          [\n"+
						"            "+l31.getLongitude()+",\n"+
						"            "+l31.getLatitude()+"\n"+
						"          ]\n"+
						"        ]\n"+
						"      }\n"+
						"    }";
		sects.add(GeoJSONFaultSection.fromFeature(Feature.fromJSON(sect3JSON)));

		double minSupraMag = 6.05d;
		double maxSingleFaultMag = 6.65d;
		//	double[] maxMultiFaultMags = { 7.95d, 7.25d };
		//	double[] maxMultiFaultMags = { 7.95d, 7.95d };
		double[] maxMultiFaultMags = { 7.15d, 7.45d };

		Range xRange = new Range(6d, 7.5d);
		Range yRange = new Range(1e-6, 1e-2);

		double maxMag = StatUtils.max(maxMultiFaultMags);
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(maxMag);
		
		ScalingRelationships scale = ScalingRelationships.HANKS_BAKUN_08;

		List<ClusterRupture> rups = new ArrayList<>();
		// add single fault ruptures
		List<Double> rupMags = new ArrayList<>();
		List<Double> rupAreas = new ArrayList<>();
		List<FaultSubsectionCluster> sectClusters = new ArrayList<>();
		for (FaultSection sect : sects)
			sectClusters.add(new FaultSubsectionCluster(List.of(sect)));
		for (FaultSubsectionCluster cluster : sectClusters) {
			for (int i=refMFD.getClosestXIndex(minSupraMag); i<=refMFD.getClosestXIndex(maxSingleFaultMag); i++) {
				double mag = refMFD.getX(i);
				double area = scale.getArea(mag, cluster.startSect.getOrigDownDipWidth()*1e3);
				ClusterRupture rup = new ClusterRupture(cluster);
				rups.add(rup);
				rupMags.add(mag);
				rupAreas.add(area);
			}
		}
		FaultSystemRupSet singleFaultRupSet = FaultSystemRupSet.builderForClusterRups(sects, new ArrayList<>(rups))
				.rupMags(Doubles.toArray(rupMags))
				.rupAreas(Doubles.toArray(rupAreas))
				.build();
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(sects);
		// now add multi-fault ruptures, sect1 -> sect2 -> sect3
		ClusterRupture baseRup = new ClusterRupture(sectClusters.get(0));
		
		FaultSubsectionCluster destCluster = sectClusters.get(1);
		double destMaxMag = maxMultiFaultMags[0];
		double dist = distAzCalc.getDistance(sects.get(0), destCluster.startSect);
		Jump jump = new Jump(sects.get(0), sectClusters.get(0), destCluster.startSect, destCluster, dist);
		ClusterRupture rup01 = baseRup.take(jump);
		for (int i=refMFD.getClosestXIndex(maxSingleFaultMag)+1; i<=refMFD.getClosestXIndex(destMaxMag); i++) {
			double mag = refMFD.getX(i);
			double area = scale.getArea(mag, destCluster.startSect.getOrigDownDipWidth()*1e3);
			rups.add(rup01);
			rupMags.add(mag);
			rupAreas.add(area);
		}
		
		FaultSystemRupSet twoFaultRupSet = FaultSystemRupSet.builderForClusterRups(sects, new ArrayList<>(rups))
				.rupMags(Doubles.toArray(rupMags))
				.rupAreas(Doubles.toArray(rupAreas))
				.build();
		
		destCluster = sectClusters.get(2);
		destMaxMag = maxMultiFaultMags[1];
		dist = distAzCalc.getDistance(sects.get(1), destCluster.startSect);
		jump = new Jump(sects.get(1), sectClusters.get(1), destCluster.startSect, destCluster, dist);
		ClusterRupture rup012 = rup01.take(jump);
		for (int i=refMFD.getClosestXIndex(maxSingleFaultMag)+1; i<=refMFD.getClosestXIndex(destMaxMag); i++) {
			double mag = refMFD.getX(i);
			double area = scale.getArea(mag, destCluster.startSect.getOrigDownDipWidth()*1e3);
			rups.add(rup012);
			rupMags.add(mag);
			rupAreas.add(area);
		}
		
		for (int r=0; r<rups.size(); r++)
			System.out.println("Rupture "+r+", M="+rupMags.get(r).floatValue()+": "+rups.get(r));

		FaultSystemRupSet fullRupSet = FaultSystemRupSet.builderForClusterRups(sects, rups)
				.rupMags(Doubles.toArray(rupMags))
				.rupAreas(Doubles.toArray(rupAreas))
				.build();
		// fake sect areas to make nucleation/participation calcs work
		// make sections have typical subsection area, half as long as wide
		double[] sectAreas = new double[sects.size()];
		for (int s=0; s<sectAreas.length; s++) {
			FaultSection sect = sects.get(s);
			double width = sect.getOrigDownDipWidth()*1e3;
			double len = width*0.5;
			sectAreas[s] = len*width;
		}
		SectAreas areas = SectAreas.precomputed(singleFaultRupSet, sectAreas);
		singleFaultRupSet.addModule(areas);
		twoFaultRupSet.addModule(areas);
		fullRupSet.addModule(areas);

		ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb.D = true;
		SectNucleationMFD_Estimator segAdjuster = new ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb(segModel, 50, true);
//		String name = "Relative G-R MFD Adjustment";
		String name = " ";
		String prefix = "multi_fault";
		
		// do this one only penalizing jump 2, as would be done when calculating it's relative GR
		double jump2Prob = segModel.calcJumpProbability(dist2);
		double jump3Prob = segModel.calcJumpProbability(dist3);
		SectNucleationMFD_Estimator penalizeTwoOnlySegAdjuster =
				new ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb(new JumpProbabilityCalc() {

					@Override
					public boolean isDirectional(boolean splayed) {
						return segModel.isDirectional(splayed);
					}

					@Override
					public String getName() {
						return segModel.getName();
					}

					@Override
					public double calcJumpProbability(ClusterRupture fullRupture, Jump jump, boolean verbose) {
						return jump2Prob;
					}
					
				}, 50, true);
		
		Color threeFaultColor = Color.RED.darker();
		Color twoFaultColor = Color.BLUE.darker();
		Color singleFaultColor = Color.CYAN.darker();
		Color fullAdjFaultColor = Color.GREEN.darker();

		System.out.println("Calculating single fault MFDs");
		IncrementalMagFreqDist singleMFD = new SupraSeisBValInversionTargetMFDs.Builder(
				singleFaultRupSet, supraB).applyDefModelUncertainties(false)
				.build().getOnFaultSupraSeisNucleationMFDs().get(0);
		
		System.out.println("Calculating two fault original MFDs");
		IncrementalMagFreqDist twoOrigMFD = new SupraSeisBValInversionTargetMFDs.Builder(
				twoFaultRupSet, supraB).applyDefModelUncertainties(false)
				.build().getOnFaultSupraSeisNucleationMFDs().get(0);
		
		System.out.println("Calculating two fault adjusted MFDs");
		IncrementalMagFreqDist penalizeTwoOnlyAdjMFD = new SupraSeisBValInversionTargetMFDs.Builder(
				fullRupSet, supraB).applyDefModelUncertainties(false)
				.adjustTargetsForData(penalizeTwoOnlySegAdjuster)
				.build().getOnFaultSupraSeisNucleationMFDs().get(0);
		
		System.out.println("Calculating three fault original MFDs");
		IncrementalMagFreqDist threeOrigMFD = new SupraSeisBValInversionTargetMFDs.Builder(
				fullRupSet, supraB).applyDefModelUncertainties(false)
				.build().getOnFaultSupraSeisNucleationMFDs().get(0);
		
		System.out.println("Calculating three fault adjusted MFDs");
		IncrementalMagFreqDist threeAdjMFD = new SupraSeisBValInversionTargetMFDs.Builder(
				fullRupSet, supraB).applyDefModelUncertainties(false)
				.adjustTargetsForData(segAdjuster)
				.build().getOnFaultSupraSeisNucleationMFDs().get(0);
		
		// figure out weights
		double w1 = 1d;
		double w2 = calcAveragingWeight(threeOrigMFD, penalizeTwoOnlyAdjMFD, singleMFD);
		double w3 = calcAveragingWeight(threeOrigMFD, threeAdjMFD, singleMFD);
//		double w2 = calcAveragingWeight(twoOrigMFD, threeAdjMFD, singleMFD);

		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DecimalFormat oDF = new DecimalFormat("0.0");
		
		singleMFD.setName("Single-Fault");
		funcs.add(singleMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, singleFaultColor));

		twoOrigMFD.setName("Two-Fault G-R");
		funcs.add(twoOrigMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, twoFaultColor));

		threeOrigMFD.setName("Three-Fault G-R");
		funcs.add(threeOrigMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, threeFaultColor));
		
		threeAdjMFD.setName("Adjusted Target");
		funcs.add(threeAdjMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, fullAdjFaultColor));
		
		boolean test = false;
		if (test) {
			IncrementalMagFreqDist testMFD = new IncrementalMagFreqDist(
					threeAdjMFD.getMinX(), threeAdjMFD.size(), threeAdjMFD.getDelta());
			
			double wDelta = w1-w2;
			for (int i=0; i<singleMFD.size(); i++)
				testMFD.add(i, singleMFD.getY(i)*wDelta);
			
			wDelta = w2-w3;
			for (int i=0; i<twoOrigMFD.size(); i++)
				testMFD.add(i, twoOrigMFD.getY(i)*wDelta);
			
			wDelta = w3;
			for (int i=0; i<threeOrigMFD.size(); i++)
				testMFD.add(i, threeOrigMFD.getY(i)*wDelta);
			
			testMFD.setName("Test");
			funcs.add(testMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, new Color(0, 0, 0, 127)));
		}
		
		addMFDOverlay(twoOrigMFD, twoFaultColor, threeAdjMFD, funcs, chars);
		addMFDOverlay(threeOrigMFD, threeFaultColor, threeAdjMFD, funcs, chars);
		
		List<IncrementalMagFreqDist> annMFDs = new ArrayList<>();
		List<Color> annColors = new ArrayList<>();
		List<String> annLabels = new ArrayList<>();
		
		annMFDs.add(singleMFD);
		annColors.add(singleFaultColor);
		annLabels.add("Single-Fault, φ₀="+oDF.format(w1));
		
		annMFDs.add(twoOrigMFD);
		annColors.add(twoFaultColor);
		annLabels.add("Two-Fault, F≤"+oDF.format(jump2Prob)+", φ₁="+oDF.format(w2));
		
		annMFDs.add(threeAdjMFD);
		annColors.add(threeFaultColor);
		annLabels.add("Three-Fault, F≤"+oDF.format(jump3Prob)+", φ₂="+oDF.format(w3));
		
		double yDeltaFract = 0.075;
		double logMaxY = Math.log10(yRange.getUpperBound());
		double logMinY = Math.log10(yRange.getLowerBound());
		double logYspan = logMaxY - logMinY;
		
		double yCount = 0.75; // start a bit down
		double yDeltaEach = logYspan*yDeltaFract;
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		for (int m=0; m<annMFDs.size(); m++) {
			IncrementalMagFreqDist annMFD = annMFDs.get(m);
			double logY = logMaxY - yDeltaEach*yCount;
			double y = Math.pow(10, logY);
			yCount++;
			
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			
			int leftBin = -1;
			int rightBin = 0;
			for (int i=0; i<annMFD.size(); i++) {
				if (annMFD.getY(i) > 0d) {
					if (leftBin < 0)
						leftBin = i;
					rightBin = i;
				}
			}
			double leftX = annMFD.getX(leftBin) - 0.4*annMFD.getDelta();
			double rightX = annMFD.getX(rightBin) + 0.4*annMFD.getDelta();
			
			double textY = Math.pow(10, Math.log10(y)+0.1);
			double topY = y;
			double botY = Math.pow(10, Math.log10(y)-0.1);
			
			xy.set(leftX, botY);
			xy.set(leftX, topY);
			xy.set(rightX, topY);
			xy.set(rightX, botY);
			
			Color annColor = annColors.get(m);
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, annColor));
			
			double textX = 0.5*(leftX+rightX);
			String label = annLabels.get(m);
			XYTextAnnotation ann = new XYTextAnnotation(label, textX, textY);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			ann.setPaint(annColor);
			ann.setTextAnchor(TextAnchor.CENTER);
			anns.add(ann);
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, name,
				"Magnitude", "Target Nucleation Rate");
		spec.setLegendVisible(true);
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);

		PlotUtils.writePlots(outputDir, prefix, gp, 800, 650, true, true, false);
	}
	
	private static void addMFDOverlay(IncrementalMagFreqDist below, Color color, IncrementalMagFreqDist above,
			List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars) {
//		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 160);
		for (int i=0; i<below.size(); i++) {
			double origY = below.getY(i);
			double adjY = above.getY(i);

			if ((float)adjY >= (float)origY) {
				// plot it on top
				double x = below.getX(i);
				double x1 = x - 0.4*below.getDelta();
				double x2 = x + 0.4*below.getDelta();
				ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
				func.set(x1, origY);
				func.set(x2, origY);
				funcs.add(func);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
			}
		}
	}
	
	private static double calcAveragingWeight(IncrementalMagFreqDist mfd, IncrementalMagFreqDist finalMFD,
			IncrementalMagFreqDist singleMFD) {
		// find a magnitude bin that exists in the given MFD but not in the single-fault MFD
		int magIndex = -1;
		for (int i=0; i<mfd.size(); i++) {
			if (mfd.getY(i) > 0d && (i >= singleMFD.size() || singleMFD.getY(i) == 0d)) {
				magIndex = i;
			}
		}
		
		System.out.println("Calculating weight based on M"+(float)mfd.getX(magIndex));
		
		double myY = mfd.getY(magIndex);
		double finalY = finalMFD.getY(magIndex);
		
		System.out.println("\tAdj Y="+(float)finalY);
		System.out.println("\tOrig Y="+(float)myY);
		
		Preconditions.checkState(finalY <= myY);
		
		System.out.println("\tWeight="+(float)finalY+" / "+(float)myY+" = "+(float)(finalY/myY));
		
		return finalY/myY;
	}

}
