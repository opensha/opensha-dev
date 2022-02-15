package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectAreas;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SegmentationImpliedSectNuclMFD_Estimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SegmentationImpliedSectNuclMFD_Estimator.MultiBinDistributionMethod;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;

public class SegAdjustmentTests {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/seg_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		// create fake rupture set
		
		double supraB = 0.8;
		Shaw07JumpDistProb segModel = new Shaw07JumpDistProb(1, 3);
		
		List<FaultSection> sects = new ArrayList<>();
		double len1 = 50d;
		double len2 = 50d;
		double len3 = 30d;
		
		double az1 = 0.5*Math.PI;
		double az2 = 0.5*Math.PI;
		double az3 = 0.25*Math.PI;
		
//		double dist2 = 10d;
//		double dist3 = 5d;
		double dist2 = segModel.calcJumpDistance(0.01);
		double dist3 = segModel.calcJumpDistance(0.3);
		
		Location l10 = new Location(34, -118);
		Location l11 = LocationUtils.location(l10, az1, len1);
		
		Location l20 = LocationUtils.location(l11, az2, dist2);
		Location l21 = LocationUtils.location(l20, az2, len2);
		
		Location l30 = LocationUtils.location(l11, az3, dist3);
		Location l31 = LocationUtils.location(l30, az3, len3);
		
//		Location l3 = new Location
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
				"        \"SlipRate\": 1\n"+
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
//		System.out.println(sect1JSON);
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
				"        \"SlipRate\": 1\n"+
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
				"        \"SlipRate\": 1\n"+
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
		double maxSingleFaultMag = 6.95d;
		double[] maxMultiFaultMags = { 7.95d, 7.25d };
		
		double maxMag = StatUtils.max(maxMultiFaultMags);
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(maxMag);
		
		int maxSingleIndex = refMFD.getClosestXIndex(maxSingleFaultMag);
		int maxMagIndex = refMFD.getClosestXIndex(maxMag);
		
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
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(sects);
		// now add multi-fault ruptures, all involving sect1
		ClusterRupture baseRup = new ClusterRupture(sectClusters.get(0));
		for (int s=0; s<sectClusters.size()-1; s++) {
			FaultSubsectionCluster destCluster = sectClusters.get(s+1);
			double destMaxMag = maxMultiFaultMags[s];
			double dist = distAzCalc.getDistance(sects.get(0), destCluster.startSect);
			Jump jump = new Jump(sects.get(0), sectClusters.get(0), destCluster.startSect, destCluster, dist);
			System.out.println("Distance from 0 to "+(s+1)+": "+dist);
			System.out.println("\t"+segModel.getName()+": "+segModel.calcJumpProbability(dist));
			for (int i=refMFD.getClosestXIndex(maxSingleFaultMag)+1; i<=refMFD.getClosestXIndex(destMaxMag); i++) {
				double mag = refMFD.getX(i);
				double area = scale.getArea(mag, destCluster.startSect.getOrigDownDipWidth()*1e3);
				rups.add(baseRup.take(jump));
				rupMags.add(mag);
				rupAreas.add(area);
			}
		}
		
		for (int r=0; r<rups.size(); r++)
			System.out.println("Rupture "+r+", M="+rupMags.get(r).floatValue()+": "+rups.get(r));
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.builderForClusterRups(sects, rups)
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
		SectAreas areas = SectAreas.precomputed(rupSet, sectAreas);
		rupSet.addModule(areas);
		
		SegmentationImpliedSectNuclMFD_Estimator segAdjuster = new SegmentationImpliedSectNuclMFD_Estimator(segModel);
		
		SupraSeisBValInversionTargetMFDs.Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, supraB);
		SupraSeisBValInversionTargetMFDs noAdjMFDs = mfdBuilder.build();
		
		mfdBuilder.adjustTargetsForData(segAdjuster);
		
		MultiBinDistributionMethod[] methods = {
				MultiBinDistributionMethod.GREEDY,
				MultiBinDistributionMethod.FULLY_DISTRIBUTED,
				MultiBinDistributionMethod.CAPPED_DISTRIBUTED
		};
		Color[] colors = {
				Color.GREEN.darker(),
				Color.BLUE,
				Color.RED.darker()
		};
		
		for (int m=0; m<methods.length; m++) {
			MultiBinDistributionMethod method = methods[m];
			Color color = colors[m];
			
			segAdjuster.setTrackIndepJumpTargets(true);
			segAdjuster.setBinDistMethod(method);
			SupraSeisBValInversionTargetMFDs adjMFDs = mfdBuilder.build();
			List<List<IncrementalMagFreqDist>> faultTargetMFDs = segAdjuster.getIndepJumpTargetSupraSeisMFDs();
			segAdjuster.setTrackIndepJumpTargets(false);
			
			for (int s=0; s<sects.size(); s++) {
				IncrementalMagFreqDist origMFD = noAdjMFDs.getOnFaultSupraSeisNucleationMFDs().get(s);
				IncrementalMagFreqDist adjMFD = adjMFDs.getOnFaultSupraSeisNucleationMFDs().get(s);
				
				List<IncrementalMagFreqDist> faultTargets = faultTargetMFDs.get(s);
				SummedMagFreqDist excess = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				for (IncrementalMagFreqDist target : faultTargets)
					excess.addIncrementalMagFreqDist(target);
				
				boolean hasExcess = false;
				for (int i=0; i<excess.size(); i++)
					if ((float)excess.getY(i) > (float)adjMFD.getY(i))
						hasExcess = true;
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				IncrementalMagFreqDist origSingleFault = origMFD.deepClone();
				for (int i=maxSingleIndex+1; i<refMFD.size(); i++)
					origSingleFault.set(i, 0d);
				IncrementalMagFreqDist origMultiFault = origMFD.deepClone();
				for (int i=0; i<=maxSingleIndex; i++)
					origMultiFault.set(i, 0d);
				
//				origMultiFault.setName("Original Multi-Fault");
//				funcs.add(origMultiFault);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, Color.LIGHT_GRAY));


				
				if (hasExcess) {
					excess.setName("Excess");
					funcs.add(excess);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, new Color(0, 255, 255, 127)));
//							new Color(color.getRed(), color.getGreen(), color.getBlue(), 127)));
				}
				
//				if (hasExcess && noRedistMFDs != null) {
//					adjMFD.setName("Redistribution");
//					funcs.add(adjMFD);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, new Color(255, 0, 255, 127)));
////							new Color(255, 0, 0, 127)));
//					
//					IncrementalMagFreqDist noRedist = noRedistMFDs.getSectSupraSeisNuclMFDs().get(s);
//					noRedist.setName("Without Redistribution");
//					funcs.add(noRedist);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, color));
//				} else {
					adjMFD.setName("Adjusted");
					funcs.add(adjMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, color));
//				}

//				origSingleFault.setName("Original Singal-Fault");
//				funcs.add(origSingleFault);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, Color.GRAY));
				
				origMFD.setName("Original");
				funcs.add(origMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, new Color(0, 0, 0, 60)));
				
				boolean first = true;
				for (EvenlyDiscretizedFunc faultTarget : faultTargets) {
					EvenlyDiscretizedFunc copy = new EvenlyDiscretizedFunc(
							faultTarget.getMinX(), faultTarget.getMaxX(), faultTarget.size());
					for (int i=0; i<faultTarget.size(); i++) {
						if (faultTarget.getY(i) == 0d)
							copy.set(i, Double.NaN);
						else
							copy.set(i, faultTarget.getY(i));
					}
					faultTarget = copy;
					if (first)
						faultTarget.setName("Fault Target"+(faultTargets.size() > 1 ? "s" : ""));
					else
						faultTarget.setName(null);
					
					funcs.add(faultTarget);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, PlotSymbol.FILLED_CIRCLE, 3f, Color.BLACK));
					first = false;
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, method.toString(),
						"Magnitude", "Target Nucleation Rate");
				spec.setLegendVisible(true);
				
				Range xRange = new Range(minSupraMag-0.05, maxMag+0.05);
				Range yRange = new Range(1e-9, 1e-4);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, false, true, xRange, yRange);
				
				PlotUtils.writePlots(outputDir, "sect_"+s+"_target_"+method.name(), gp, 800, 650, true, false, false);
			}
		}
	}

}
