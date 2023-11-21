package scratch.kevin.nshm23.ceusFSS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.ned.nshm23.CEUS_FSS_creator;
import scratch.ned.nshm23.CEUS_FSS_creator.FaultModelEnum;

public class NedConverterMerge {

	public static void main(String[] args) throws IOException {
		String path = "/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0/";
		
		File outputDir = new File("/data/kevin/nshm23/ceus_converted_fss");
		
		for (FaultModelEnum fm : FaultModelEnum.values()) {
			System.out.println("Building for "+fm);
			List<FaultSystemSolution> sols = CEUS_FSS_creator.getFaultSystemSolutionList(path, fm);
			System.out.println("Fetched "+sols.size()+" solutions");
			
			FaultSystemSolution merged = null;
			if (sols.size() == 1) {
				merged = sols.get(0);
				merged.write(new File(outputDir, "sol_"+fm.name()+".zip"));
			} else {
				// merge them
				int totNumSects = 0;
				int totNumRups = 0;
				for (int i=0; i<sols.size(); i++) {
					FaultSystemSolution sol = sols.get(i);
					FaultSystemRupSet rupSet = sol.getRupSet();
					System.out.println("RupSet "+i+" has "+rupSet.getNumSections()+" sects, "+rupSet.getNumRuptures()+" rups");
					totNumSects += rupSet.getNumSections();
					totNumRups += rupSet.getNumRuptures();
					sol.write(new File(outputDir, "sol_"+fm.name()+"_"+i+".zip"));
				}
				System.out.println("Total: "+totNumSects+" sects, "+totNumRups+" rups");
				
				List<FaultSection> mergedSects = new ArrayList<>(totNumSects);
				List<List<Integer>> sectionForRups = new ArrayList<>(totNumSects);
				double[] mags = new double[totNumRups];
				double[] rakes = new double[totNumRups];
				double[] rupAreas = new double[totNumRups];
				double[] rupLengths = new double[totNumRups];
				double[] rates = new double[totNumRups];
				
				int sectIndex = 0;
				int rupIndex = 0;
				
				for (FaultSystemSolution sol : sols) {
					FaultSystemRupSet rupSet = sol.getRupSet();
					int[] sectMappings = new int[rupSet.getNumSections()];
					System.out.println("Merging sol with "+rupSet.getNumSections()+" sects and "+rupSet.getNumRuptures()+" rups");
					for (int s=0; s<sectMappings.length; s++) {
						FaultSection sect = rupSet.getFaultSectionData(s);
						sect = sect.clone();
						sectMappings[s] = sectIndex;
						sect.setSectionId(sectIndex);
						mergedSects.add(sect);
						
						sectIndex++;
					}
					
					for (int r=0; r<rupSet.getNumRuptures(); r++) {
						List<Integer> prevSectIDs = rupSet.getSectionsIndicesForRup(r);
						List<Integer> newSectIDs = new ArrayList<>(prevSectIDs.size());
						for (int s : prevSectIDs)
							newSectIDs.add(sectMappings[s]);
						sectionForRups.add(newSectIDs);
						mags[rupIndex] = rupSet.getMagForRup(r);
						rakes[rupIndex] = rupSet.getAveRakeForRup(r);
						rupAreas[rupIndex] = rupSet.getAreaForRup(r);
						rupAreas[rupIndex] = rupSet.getLengthForRup(r);
						rates[rupIndex] = sol.getRateForRup(r);
						
						rupIndex++;
					}
				}
				
				FaultSystemRupSet mergedRupSet = new FaultSystemRupSet(mergedSects, sectionForRups, mags, rakes, rupAreas, rupLengths);
				merged = new FaultSystemSolution(mergedRupSet, rates);
				
				merged.write(new File(outputDir, "sol_"+fm.name()+"_merged.zip"));
			}
		}
	}

}
