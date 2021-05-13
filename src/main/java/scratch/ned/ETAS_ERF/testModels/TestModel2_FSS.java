/**
 * 
 */
package scratch.ned.ETAS_ERF.testModels;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

/**
 * @author field
 *
 */
public class TestModel2_FSS extends FaultSystemSolution {
	
	final static boolean D = false;	// debug flag
	
	int minNumSectInRup=4;
	double rake=0;
	double slipRate = 25;	// mm/yr
	double ddw=12;
	double dip=90;
	double minLonForFault=241;
//	double maxLonForFault=243;
	double maxLonForFault=249;
	Location faultEndLoc1 = new Location(36,minLonForFault-360);
	Location faultEndLoc2 = new Location(36,maxLonForFault-360);
	
	int totNumRups;	// found by computing once
	
	ArrayList<FaultSectionPrefData> subSectionData;
	double[] rateForRup;
	double[] magForRup;
	double[] areaForRup;	// square-meters
	double[] aveSlipForRup;	// meters
	double[] distFromCenterForRupWeight;	// meters
	int[] firstSubSectForRup;
	int[] lastSubSectForRup;
	
	
	// MFDs
	GutenbergRichterMagFreqDist targetFaultGR;	// goes to M 2.5
	ArbIncrementalMagFreqDist faultGR;			// goes to the smallest fault mag
	
	
	public TestModel2_FSS() {
		super();
		
		FaultTrace trace = new FaultTrace(null);
		trace.add(faultEndLoc1);
		trace.add(faultEndLoc2);
		FaultSectionPrefData faultSectData = new FaultSectionPrefData();
		faultSectData.setAseismicSlipFactor(0.0);
		faultSectData.setAveDip(dip);
		faultSectData.setAveLowerDepth(ddw);
		faultSectData.setAveUpperDepth(0);
		faultSectData.setAveSlipRate(slipRate);
		faultSectData.setFaultTrace(trace);
		
		double subSectionMaxLength = ddw/minNumSectInRup;
		
		double area = faultSectData.getTraceLength()*faultSectData.getReducedDownDipWidth();  // sq-km
		
		subSectionData = faultSectData.getSubSectionsList(subSectionMaxLength);
		
		double[] areaForSections = new double[subSectionData.size()];	// square-meters
		for(int s=0;s<subSectionData.size();s++) {
			subSectionData.get(s).setSectionName("Subsection "+s);
			subSectionData.get(s).setParentSectionName("Bogus");
			subSectionData.get(s).setSectionId(s);
//			subSectionData.get(s).setDateOfLastEvent(0);
			areaForSections[s] = subSectionData.get(s).getTraceLength()*subSectionData.get(s).getReducedDownDipWidth();
//			if(D)
//				System.out.println("subsection name = "+subSectionData.get(s).getName());
		}

		
		FaultSectionPrefData firstSubSec = subSectionData.get(0);
		
		
		if(D) {
			System.out.println("subsection lengths = "+(float)firstSubSec.getTraceLength()+
					"; num subSections = "+subSectionData.size());			
		}
		
		
//		Ellsworth_B_WG02_MagAreaRel ellB_magArea = new Ellsworth_B_WG02_MagAreaRel();
		HanksBakun2002_MagAreaRel hbMagArea = new HanksBakun2002_MagAreaRel();
		double maxMag = hbMagArea.getMedianMag(area);
		double minMag = hbMagArea.getMedianMag(firstSubSec.getReducedDownDipWidth()*minNumSectInRup*firstSubSec.getTraceLength());

		if(D) {
			System.out.println("\nminMag="+(float)minMag+"; maxMag="+(float)maxMag);
		}
		
		double totMoRate = FaultMomentCalc.getMoment(area*1e6, slipRate*1e-3);
		
		totMoRate *= 0.5/1.1420689;	// this makes it one M>=5 every 2 years
//		totMoRate /= 10;
		
		EvenlyDiscretizedFunc magNumDist = new EvenlyDiscretizedFunc(5.15, 9.55, 45);
		magNumDist.setTolerance(0.1);
		
		// compute tot num ruptures
		totNumRups=0;
		for(int s=minNumSectInRup;s<subSectionData.size()+1;s++) {
			int numSectForRup=s;
			if(numSectForRup<subSectionData.size()/2)
				totNumRups += subSectionData.size()+1 - s;
		}
		rateForRup = new double[totNumRups];
		magForRup = new double[totNumRups];
		areaForRup = new double[totNumRups];	// square-meters
		aveSlipForRup = new double[totNumRups];	// meters
		firstSubSectForRup = new int[totNumRups];
		lastSubSectForRup = new int[totNumRups];
		distFromCenterForRupWeight = new double[totNumRups];
		
		int rupIndex=0;
		for(int s=minNumSectInRup;s<subSectionData.size()+1;s++) {
			int numSectForRup=s;
			if(numSectForRup<subSectionData.size()/2) {
				int numRups = subSectionData.size()+1 - s;
				double length = (s)*firstSubSec.getTraceLength();
				double rupArea = length*firstSubSec.getReducedDownDipWidth();
				double mag = hbMagArea.getMedianMag(rupArea);
				mag = ((double)Math.round(mag*100))/100;
				double aveSlip = FaultMomentCalc.getSlip(rupArea*1e6, MagUtils.magToMoment(mag));
				magNumDist.add(mag, 1.0);	// use this if weighting by distance from end of fault
//				magNumDist.add(mag, numRups);
				if(D) {
					System.out.println("\tMag="+(float)mag+" for "+s+" sub sections; numRups="+numRups);
				}
				int tempRupIndex = rupIndex;
				double totWeight=0;
				for(int r=0;r<numRups;r++) {
					magForRup[rupIndex] = mag;
					aveSlipForRup[rupIndex]=aveSlip;
					areaForRup[rupIndex]=rupArea*1e6;	// converted to meters-squared
					firstSubSectForRup[rupIndex]=r;
					lastSubSectForRup[rupIndex]=r+s-1;
					double relDistFromCenter = Math.abs(((double)lastSubSectForRup[rupIndex]+(double)firstSubSectForRup[rupIndex])/2d - (double)subSectionData.size()/2d);	// rup midpoint minus fault center
					relDistFromCenter = (relDistFromCenter)/((double)subSectionData.size()/2d);
					if(relDistFromCenter>=1.0)
						throw new RuntimeException();
					if(firstSubSectForRup[rupIndex]==0 || lastSubSectForRup[rupIndex] == subSectionData.size()-1)
						distFromCenterForRupWeight[rupIndex] = 1e-2;
					else
						distFromCenterForRupWeight[rupIndex] = 1.0;

//					distFromCenterForRupWeight[rupIndex] = Math.sqrt(Math.cos(relDistFromCenter*Math.PI/2d));
//					distFromCenterForRupWeight[rupIndex] = Math.pow(10,-Math.abs(relDistFromCenter)/30d);	// ~10-4 reduction at edges
					totWeight += distFromCenterForRupWeight[rupIndex];
					if(D) {
						System.out.println("\t\t"+(float)magForRup[rupIndex]+"\t"+Math.round(areaForRup[rupIndex])+"\t"+
								(float)aveSlipForRup[rupIndex]+"\t"+firstSubSectForRup[rupIndex]+"\t"+lastSubSectForRup[rupIndex]+"\t"+relDistFromCenter);
					}
					rupIndex+=1;
				}
				// normalize weights so they sum to 1.0
				double testTotWt=0;
				for(int r=0;r<numRups;r++) {
					distFromCenterForRupWeight[tempRupIndex+r] /= totWeight;
					testTotWt+=distFromCenterForRupWeight[tempRupIndex+r];
				}
				if(D) {
					System.out.println("\t\ttestTotWt="+(float)testTotWt);
				}
			}
		}
		
		
		if(D) {
			System.out.println("\ntotNumRups="+totNumRups);
			System.out.println("\nmagForRup.length="+magForRup.length);
			System.out.println("\n"+magNumDist);
		}
		
		double roundedMaxMag = magNumDist.getX(magNumDist.getClosestXIndex(maxMag));
		double roundedMinMag = magNumDist.getX(magNumDist.getClosestXIndex(minMag));
		int numPts1 = (int)Math.round((roundedMaxMag-2.55)/0.1) + 1;
		
		targetFaultGR = new GutenbergRichterMagFreqDist(2.55, numPts1, 0.1,2.55, roundedMaxMag,totMoRate, 1.0);
		targetFaultGR.setName("targetFaultGR");
//		System.out.println("\n"+targetFaultGR);
		
		int numPts2 = (int)Math.round((roundedMaxMag-roundedMinMag)/0.1) + 1;

		faultGR = new ArbIncrementalMagFreqDist(roundedMinMag, roundedMaxMag, numPts2);
		faultGR.setName("faultGR");
		for(int i=0;i<faultGR.size();i++) {
			double mag = faultGR.getX(i);
			faultGR.set(mag, targetFaultGR.getY(mag));
		}
		
		ArbIncrementalMagFreqDist testMFD = new ArbIncrementalMagFreqDist(roundedMinMag, roundedMaxMag, numPts2);
		testMFD.setName("testMFD");
		testMFD.setInfo("Implied by model; should equal faultGR");

		for(int r=0; r<totNumRups;r++) {
			int index = faultGR.getClosestXIndex(magForRup[r]);
			rateForRup[r] = distFromCenterForRupWeight[r]*faultGR.getY(index)/magNumDist.getY(magForRup[r]);	// apply dist from center weight
//			rateForRup[r] = faultGR.getY(index)/magNumDist.getY(magForRup[r]);
			testMFD.add(index, rateForRup[r]);
		}

		
//		System.out.println("\n"+faultGR);
		
		if(D) {
			ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
			funcs.add(faultGR);
			funcs.add(faultGR.getCumRateDistWithOffset());
			funcs.add(targetFaultGR);
			funcs.add(targetFaultGR.getCumRateDistWithOffset());
			funcs.add(testMFD);

			GraphWindow graph = new GraphWindow(funcs, ""); 

			System.out.println("MomentRates: "+totMoRate+"\t"+faultGR.getTotalMomentRate()
					+"\t"+targetFaultGR.getTotalMomentRate());

		}
		
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		for (int r=0; r<totNumRups; r++) {
			ArrayList<Integer> indices = new ArrayList<Integer>();
			for(int i=firstSubSectForRup[r]; i<=lastSubSectForRup[r]; i++)
				indices.add(i);
			sectionForRups.add(indices);
		}
		
		double[] rakes = new double[totNumRups];
		
		FaultSystemRupSet rupSet = new FaultSystemRupSet(subSectionData, null, null, areaForSections, sectionForRups, magForRup,
				rakes, areaForRup, null, null);
		
		System.out.println("rupSet.getNumRuptures()="+rupSet.getNumRuptures());
		
		init(rupSet, rateForRup, null, null);
		
		if(D) {
			double[] partRateArray = calcParticRateForAllSects(0d, 10d);
			EvenlyDiscretizedFunc partRateFunc = new EvenlyDiscretizedFunc(0d,partRateArray.length,1.0);
			for(int i=0;i<partRateArray.length;i++)
				partRateFunc.set(i,partRateArray[i]);
			partRateFunc.setName("Section Participation Rates");
			partRateFunc.setInfo("Max/Min="+(float)(partRateFunc.getMaxY()/partRateFunc.getMinY()));
			GraphWindow graph2 = new GraphWindow(partRateFunc, "Section Participation Rates");
		}
	}
	
	/**
	 * This is the target MFD for the fault (going down to M 2.5)
	 * @return
	 */
	public GutenbergRichterMagFreqDist getTargetFaultGR() {
		return targetFaultGR;
	}
	
	public ArbIncrementalMagFreqDist getFaultGR() {
		return faultGR;
	}
	

	/* (non-Javadoc)
	 * @see scratch.UCERF3.FaultSystemSolution#getRateForAllRups()
	 */
	@Override
	public double[] getRateForAllRups() {
		return rateForRup;
	}

	/* (non-Javadoc)
	 * @see scratch.UCERF3.FaultSystemSolution#getRateForRup(int)
	 */
	@Override
	public double getRateForRup(int rupIndex) {
		return rateForRup[rupIndex];
	}
	
	/**
	 * This returns all ruptures that are: 1) completely inside the given rupture surface, 
	 * 2) completely surrounding it; 3) have half or more of its surface within, or 
	 * 4) extends inside by at least the given numSectOverlap.
	 * @param rthRup
	 * @param numSectOverlap
	 * @return
	 */
	public List<Integer> getRupsThatOverlapGivenRup(int rthRup, int numSectOverlap) {
		
		ArrayList<Integer> rupsWithOverlap = new ArrayList<Integer> ();
		
		List<Integer>  rupSects = getRupSet().getSectionsIndicesForRup(rthRup);

		int firstRupSect = rupSects.get(0);
		int lastRupSect = rupSects.get(rupSects.size()-1);
// System.out.println("TARGET: "+rthRup+"\t"+firstRupSect+"\t"+lastRupSect);

				
		for(int i=0; i<getRupSet().getNumRuptures(); i++) {
			List<Integer>  sects = getRupSet().getSectionsIndicesForRup(i);
			int first = sects.get(0);
			int last = sects.get(sects.size()-1);
			
			// check if it's entirely inside
			if(first >= firstRupSect && last <= lastRupSect) {
				rupsWithOverlap.add(i);
			} 
			// check if surrounding
			else if(first <= firstRupSect && last >= lastRupSect) {
				rupsWithOverlap.add(i);
			}
			
			// now examine fractional overlaps
			else {
				int numInside =0;
				for(Integer s:sects) {
					if(s>=firstRupSect && s<=lastRupSect)
						numInside +=1;
				}
				double fractInside = (double)numInside/(double)sects.size();
				if(fractInside >=0.49999999)	// if more than half the rupture is inside
					rupsWithOverlap.add(i);
				else if(numInside>=numSectOverlap)	// if more than numSectOverlap are inside
					rupsWithOverlap.add(i);
			}
			
			// write test
// System.out.println(i+"\t"+first+"\t"+last+"\t"+rupsWithOverlap.contains(i));

		}
		

		return rupsWithOverlap;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		TestModel2_FSS test = new TestModel2_FSS();
		
		System.out.println(test.getRupSet().getNumRuptures());

//		test.getRupsThatOverlapGivenRup(892, 10);
		
//		File file = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR,"/TestModel1_FSS.zip");
//		
//		try {
//			FaultSystemIO.writeSol(test, file);
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
//		for(int i=0; i<test.getNumRuptures(); i++) {
//			List<Integer>  sects = test.getSectionsIndicesForRup(i);
//			System.out.println(i+"\t"+sects.get(0)+"\t"+sects.get(sects.size()-1));
//		}

	}

}
