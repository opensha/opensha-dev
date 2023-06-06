package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
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
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;

public class SlipRateSegAdjTests {
	
	public static void main(String[] args) throws IOException {
		double supraB = 0.5;
		Shaw07JumpDistProb segModel = Shaw07JumpDistProb.forHorzOffset(1d, 3, 2);
		
		List<FaultSection> sects = new ArrayList<>();
		double len1 = 50d;
		double len2 = 50d;
		double len3 = 50d;
		
		double az1 = 0d;
		double az2 = 0d;
		double az3 = 0d;
		
//		double dist2 = 10d;
//		double dist3 = 5d;
		double dist2 = 0d;
//		double dist2 = segModel.calcJumpDistance(0.1);
		double dist3 = segModel.calcJumpDistance(0.1);
		
		double slip1 = 10d;
		double slip2 = 10d;
		double slip3 = 100d;
		
		Location l10 = new Location(34, -118);
		Location l11 = LocationUtils.location(l10, az1, len1);
		
		Location l20 = LocationUtils.location(l11, az2, dist2);
		Location l21 = LocationUtils.location(l20, az2, len2);
		
		Location l30 = LocationUtils.location(l21, az3, dist3);
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
				"        \"SlipRate\": "+slip1+"\n"+
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
				"        \"SlipRate\": "+slip2+"\n"+
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
				"        \"SlipRate\": "+slip3+"\n"+
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
//		double[] maxMultiFaultMags = { 7.95d, 7.25d };
//		double[] maxMultiFaultMags = { 7.95d, 7.95d };
		double[] maxMultiFaultMags = { 7.45d, 7.95d };
		
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
		singleFaultRupSet.addModule(areas);
		
		RelGRWorstJumpProb.D = true;
		RelGRWorstJumpProb adj = new RelGRWorstJumpProb(segModel, 50, true);
		
		SupraSeisBValInversionTargetMFDs.Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, supraB);
		mfdBuilder.applyDefModelUncertainties(false);
		
		mfdBuilder.clearTargetAdjustments().adjustTargetsForData(adj);
		
		SupraSeisBValInversionTargetMFDs mfds = mfdBuilder.build();
		
		IncrementalMagFreqDist mfd0 = mfds.getOnFaultSupraSeisNucleationMFDs().get(0);
		
		System.out.println("MFD for 0");
		System.out.println(mfd0);
		
		List<IncrementalMagFreqDist> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(mfd0);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, new Range(6d, 8d), new Range(1e-7, 1e-3));
		
		PlotUtils.writePlots(new File("/tmp"), "seg_slip_tests", gp, 800, 750, true, false, false);
	}

}
