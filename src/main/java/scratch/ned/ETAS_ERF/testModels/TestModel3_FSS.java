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
 * This one is perfectly segmented along strike
 * @author field
 *
 */
public class TestModel3_FSS extends FaultSystemSolution {
	
	final static boolean D = true;	// debug flag
	
	int numIdenticalRups = 5;  // these are identical ruptures - exact same place
	
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
	
	int totNumRups;
	
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
	
	
	public TestModel3_FSS() {
		super();
		
		double aveRupRI = 100;	// years
		
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
		
		// compute tot num ruptures
		int numUniqueRups = subSectionData.size()/minNumSectInRup;
		totNumRups=numIdenticalRups*numUniqueRups;
		
		rateForRup = new double[totNumRups];
		magForRup = new double[totNumRups];
		areaForRup = new double[totNumRups];	// square-meters
		aveSlipForRup = new double[totNumRups];	// meters
		firstSubSectForRup = new int[totNumRups];
		lastSubSectForRup = new int[totNumRups];
		
		double rupArea = minNumSectInRup*firstSubSec.getTraceLength()*firstSubSec.getReducedDownDipWidth();
		double rupMag = hbMagArea.getMedianMag(rupArea);
		double aveSlip = FaultMomentCalc.getSlip(rupArea*1e6, MagUtils.magToMoment(rupMag));
		if (D) 
			System.out.println("rupMag="+rupMag+";  aveSlip="+aveSlip);
		
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		
		for(int r=0;r<numUniqueRups;r++) {
			rateForRup[r]= ((double)(1+r)/((double)totNumRups/2))/aveRupRI;	// linear trend to get variation
			magForRup[r]=rupMag;
			areaForRup[r]=rupArea;	// square-meters
			aveSlipForRup[r]=aveSlip;	// meters
			firstSubSectForRup[r] = r*minNumSectInRup;
			lastSubSectForRup[r] = r*minNumSectInRup+3;
			ArrayList<Integer> indices = new ArrayList<Integer>();
			for(int i=firstSubSectForRup[r]; i<=lastSubSectForRup[r]; i++)
				indices.add(i);
			sectionForRups.add(indices);
			if (D) System.out.println("\t"+firstSubSectForRup[r]+"\t"+lastSubSectForRup[r]);
			
		}
		
		// make duplicates
		for(int i =1;i<numIdenticalRups;i++) {
			for(int r=0;r<numUniqueRups;r++) {
				int r_dupl = r+i*numUniqueRups;
				rateForRup[r_dupl]= rateForRup[r];	// linear trend to get variation
				magForRup[r_dupl]=magForRup[r];
				areaForRup[r_dupl]=areaForRup[r];	// square-meters
				aveSlipForRup[r_dupl]=aveSlipForRup[r];	// meters
				firstSubSectForRup[r_dupl] = firstSubSectForRup[r];
				lastSubSectForRup[r_dupl] = lastSubSectForRup[r];
				sectionForRups.add(sectionForRups.get(r));
			}
		}
		
		if(D) {
			System.out.println("\ntotNumRups="+totNumRups+"; numIdenticalRups="+numIdenticalRups);
			System.out.println("\nmagForRup.length="+magForRup.length);
		}
				
		double[] rakes = new double[totNumRups];	// all zeros
		
		FaultSystemRupSet rupSet = new FaultSystemRupSet(subSectionData, null, null, areaForSections, sectionForRups, magForRup,
				rakes, areaForRup, null, null);
		
		System.out.println("rupSet.getNumRuptures()="+rupSet.getNumRuptures());
		
		init(rupSet, rateForRup, null, null);
		
		double[] partRateArray = calcParticRateForAllSects(0d, 10d);
//		for(int s=0;s<subSectionData.size();s++) {
//			subSectionData.get(s).setDateOfLastEvent(0.5/partRateArray[s]);
//		}

		if(D) {
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
		
		TestModel3_FSS test = new TestModel3_FSS();
		
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
