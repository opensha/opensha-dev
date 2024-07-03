package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import org.opensha.commons.calc.GaussianDistCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

import scratch.kevin.pointSources.PointSourceHazardComparison.PointSourceCalcERF;
import scratch.kevin.pointSources.PointSourceHazardComparison.PointSourceType;

public class HazardTargetedDistCorrCalc {

	public static void main(String[] args) throws IOException {
		// calculate the rJB for each mag/horz dist that gives you equivalent ground motion exceedance probs using
		// PointSourceNshm as you would doing it right
		
		// we'll do this for a single period
		double period = 0d;
		// and for this conditional probability of exceedance
		double condExceedProb = 0.5;
//		double condExceedProb = GaussianDistCalc.getExceedProb(1d); // +1 sigma
		// and for this GMM
		AttenRelRef gmmRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		
		Location centerLoc = new Location(0d, 0d);
		
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(5.05, 30, 0.1);
		mfd.setAllButTotMoRate(mfd.getMinX(), mfd.getMaxX(), 1d, 1d);
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, 301, 1d);
		int numEachDist = 10;
		List<List<Location>> distLocs = new ArrayList<>();
		for (int i=0; i<distFunc.size(); i++) {
			double dist = distFunc.getX(i);
			if (dist == 0d) {
				distLocs.add(List.of(centerLoc));
			} else {
				List<Location> locs = new ArrayList<>(numEachDist);
				double deltaEach = 2d*Math.PI/numEachDist;
				for (int j=0; j<numEachDist; j++)
					locs.add(LocationUtils.location(centerLoc, deltaEach*j, dist));
				distLocs.add(locs);
			}
		}
		
		PointSourceType refType = PointSourceType.OCT_QUAD;
		int numPerRefType = 1000;
		PointSourceType targetType = PointSourceType.SIMPLE_APPROX_FINITE_POINT_NO_CORR;
		int numPerTargetType = 1;
		
		File mainDir = new File("/data/kevin/markdown/nshm23-misc/point_source_corr/");
		File subDir = new File(mainDir, targetType.name()+"_vs_"+refType.name());
		Preconditions.checkState(subDir.exists() || subDir.mkdir());
		File outputDir = new File(subDir, "hazard_equiv_rJBs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Random rand = new Random(123456789l);
		
		FocalMech[] mechs = { FocalMech.STRIKE_SLIP, FocalMech.REVERSE };
		
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		
		for (FocalMech mech : mechs) {
			System.out.println("Building ERFs for "+mech);
			PointSourceCalcERF refERF = new PointSourceCalcERF(refType, centerLoc, mfd,
					mech.rake(), mech.dip(), numPerRefType, rand);
			PointSourceCalcERF targetERF = new PointSourceCalcERF(targetType, centerLoc, mfd,
					mech.rake(), mech.dip(), numPerTargetType, rand);
			
			System.out.println("Calculating rJB as a function of center location for "+targetType);
			EvenlyDiscrXYZ_DataSet targetRJBs = calcRJBs(targetERF, distFunc, distLocs, mfd);
//			printTable(targetRJBs);
			
			System.out.println("Calculating IMs at condExceedProb="+(float)condExceedProb+" for "+refType);
			EvenlyDiscrXYZ_DataSet refIMs = calcIMsAtExceedProb(refERF, distFunc, distLocs, mfd, gmmRef, period, gmmDeque, condExceedProb);
//			printTable(refIMs);
			System.out.println("Calculating IMs at condExceedProb="+(float)condExceedProb+" for "+targetType);
			EvenlyDiscrXYZ_DataSet targetIMs = calcIMsAtExceedProb(targetERF, distFunc, distLocs, mfd, gmmRef, period, gmmDeque, condExceedProb);
//			printTable(targetIMs);
			
			EvenlyDiscrXYZ_DataSet corrRJBs = new EvenlyDiscrXYZ_DataSet(distFunc.size(), mfd.size(),
					distFunc.getMinX(), mfd.getMinX(), distFunc.getDelta(), mfd.getDelta());
			for (int m=0; m<mfd.size(); m++) {
				double minRJB = targetRJBs.get(0, m);
				double maxRJB = targetRJBs.get(distFunc.size()-1, m);
				double maxIM = targetIMs.get(0, m);
				double minIM = targetIMs.get(distFunc.size()-1, m);
				EvenlyDiscretizedFunc targetIMsFunc = new EvenlyDiscretizedFunc(distFunc.getMinX(), distFunc.size(), distFunc.getDelta());
				for (int d=0; d<distFunc.size(); d++)
					targetIMsFunc.set(d, targetIMs.get(d, m));
				double mag = mfd.getX(m);
				Preconditions.checkState(maxIM > minIM, "maxIM=%s should be greater than minIM=%s; mag=%s", maxIM, minIM, mag);
				for (int d=0; d<distFunc.size(); d++) {
					double refIM = refIMs.get(d, m);
					// find the target rJB at which targetIM matches refIM
					double corrRJB;
					if (refIM > maxIM) {
						// this value is greater than any of our IMs, use our smallest rJB
						corrRJB = minRJB;
					} else if (refIM < minIM) {
						// this value is less than any of our IMs, use our largest rJB
						corrRJB = maxRJB;
					} else {
						// interpolate
						double horzDistForIM = targetIMsFunc.getFirstInterpolatedX(refIM);
						corrRJB = targetRJBs.bilinearInterpolation(horzDistForIM, mag);
					}
					corrRJBs.set(d, m, corrRJB);
				}
			}
			
			printTable(corrRJBs);
			
			File csvFile = new File(outputDir, mech.name()+"_equiv_rJBs.csv");
			System.out.println("Writing "+csvFile.getAbsolutePath());
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>(mfd.size()+1);
			header.add("Horizontal Distance");
			for (int m=0; m<mfd.size(); m++)
				header.add((float)mfd.getX(m)+"");
			csv.addLine(header);
			for (int d=0; d<distFunc.size(); d++) {
				List<String> line = new ArrayList<>(header.size());
				line.add((float)distFunc.getX(d)+"");
				for (int m=0; m<mfd.size(); m++)
					line.add((float)corrRJBs.get(d, m)+"");
				csv.addLine(line);
			}
			
			csv.writeToFile(csvFile);
		}
	}
	
	private static EvenlyDiscrXYZ_DataSet calcRJBs(PointSourceCalcERF erf, EvenlyDiscretizedFunc distFunc,
			List<List<Location>> distLocs, IncrementalMagFreqDist mfd) {
		double[][] rJBs = new double[mfd.size()][distFunc.size()]; // y is outer index
		
		IntStream.range(0, distFunc.size()).parallel().forEach(d -> {
			int[] numEntries = new int[mfd.size()];
			
			for (Location loc : distLocs.get(d)) {
				for (ProbEqkSource source : erf) {
					for (ProbEqkRupture rup : source) {
						double mag = rup.getMag();
						int m = mfd.getClosestXIndex(mag);
						RuptureSurface surf = rup.getRuptureSurface();
						rJBs[m][d] += surf.getDistanceJB(loc);
						numEntries[m]++;
					}
				}
			}
			
			for (int m=0; m<mfd.size(); m++)
				rJBs[m][d] /= (double)numEntries[m];
		});
		
		return new EvenlyDiscrXYZ_DataSet(rJBs, distFunc.getMinX(), mfd.getMinX(), distFunc.getDelta(), mfd.getDelta());
	}
	
	private static EvenlyDiscrXYZ_DataSet calcIMsAtExceedProb(PointSourceCalcERF erf, EvenlyDiscretizedFunc distFunc,
			List<List<Location>> distLocs, IncrementalMagFreqDist mfd, AttenRelRef gmmRef, double period,
			Deque<ScalarIMR> gmmDeque, double condExceedProb) {
		double[][] ims = new double[mfd.size()][distFunc.size()]; // y is outer index
		
		List<Integer> doneIndexes = new ArrayList<>(distFunc.size());
		IntStream.range(0, distFunc.size()).parallel().forEach(d -> {
			int[] numEntries = new int[mfd.size()];
			
			ScalarIMR gmm = null;
			synchronized (gmmDeque) {
				if (!gmmDeque.isEmpty())
					gmm = gmmDeque.pop();
			}
			if (gmm == null)
				gmm = gmmRef.get();
			if (period == 0d) {
				gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				Preconditions.checkState(period > 0d);
				gmm.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			}
			
			for (Location loc : distLocs.get(d)) {
				Site site = new Site(loc);
				site.addParameterList(gmm.getSiteParams());
				gmm.setSite(site);
				for (ProbEqkSource source : erf) {
					for (ProbEqkRupture rup : source) {
						gmm.setEqkRupture(rup);
						double mag = rup.getMag();
						int m = mfd.getClosestXIndex(mag);
						ims[m][d] += gmm.getIML_AtExceedProb(condExceedProb);
						numEntries[m]++;
					}
				}
			}
			
			for (int m=0; m<mfd.size(); m++)
				ims[m][d] /= (double)numEntries[m];
			
			synchronized (gmmDeque) {
				gmmDeque.push(gmm);
				doneIndexes.add(d);
				System.out.print(".");
				if (doneIndexes.size() % 100 == 0)
					System.out.println(" "+doneIndexes.size());
			}
		});
		if (doneIndexes.size() % 100 != 0)
			System.out.println(" "+doneIndexes.size());
		
		return new EvenlyDiscrXYZ_DataSet(ims, distFunc.getMinX(), mfd.getMinX(), distFunc.getDelta(), mfd.getDelta());
	}
	
	private static DecimalFormat tableDF = new DecimalFormat("0.00");
	
	private static void printTable(EvenlyDiscrXYZ_DataSet table) {
		StringBuilder buf = new StringBuilder();
		buf.append("d/m");
		for (int m=0; m<table.getNumY(); m++)
			buf.append("\t").append(tableDF.format(table.getY(m)));
		System.out.println(buf);
		for (int d=0; d<table.getNumX(); d++) {
			buf = new StringBuilder();
			buf.append(tableDF.format(table.getX(d)));
			for (int m=0; m<table.getNumY(); m++)
				buf.append("\t").append(tableDF.format(table.get(d, m)));
			System.out.println(buf);
		}
	}

}
