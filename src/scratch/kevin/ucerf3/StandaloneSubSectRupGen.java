package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.SectionCluster;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.SectionConnectionStrategy;
import scratch.UCERF3.inversion.UCERF3SectionConnectionStrategy;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesRecord;
import scratch.UCERF3.inversion.laughTest.UCERF3PlausibilityConfig;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class StandaloneSubSectRupGen {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws DocumentException, IOException {
		// this is the input fault section data file
		File fsdFile = new File("dev/scratch/UCERF3/data/FaultModels/FM3_1.xml");
		// directory to write output files
		File outputDir = new File("/tmp");
		// maximum sub section length (in units of DDW)
		double maxSubSectionLength = 0.5;
		// max distance for linking multi fault ruptures, km
		double maxDistance = 5d;
		boolean coulombFilter = true;
		FaultModels fm = FaultModels.FM3_1;
		// this is a list of parent fault sections to remove. can be empty or null
		// currently set to remove Garlock to test Coulomb remapping.
		List<Integer> sectsToRemove = Lists.newArrayList(49, 341);
		// this is a list of sections to keep. if non null and non empty, only these
		// ids will be kept
		List<Integer> sectsToKeep = Lists.newArrayList(32);
		
		Preconditions.checkState(!coulombFilter || fm == FaultModels.FM3_1 || fm == FaultModels.FM3_2);
		
		// load in the fault section data ("parent sections")
		List<FaultSection> fsd = FaultModels.loadStoredFaultSections(fsdFile);
		
		if (sectsToRemove != null && !sectsToRemove.isEmpty()) {
			System.out.println("Removing these parent fault sections: "
					+Joiner.on(",").join(sectsToRemove));
			// iterate backwards as we will be removing from the list
			for (int i=fsd.size(); --i>=0;)
				if (sectsToRemove.contains(fsd.get(i).getSectionId()))
					fsd.remove(i);
		}
		
		if (sectsToKeep != null && !sectsToKeep.isEmpty()) {
			System.out.println("Only keeping these parent fault sections: "
					+Joiner.on(",").join(sectsToKeep));
			// iterate backwards as we will be removing from the list
			for (int i=fsd.size(); --i>=0;)
				if (!sectsToKeep.contains(fsd.get(i).getSectionId()))
					fsd.remove(i);
		}
		
		System.out.println(fsd.size()+" Parent Fault Sections");
		
		// this list will store our subsections
		List<FaultSection> subSections = Lists.newArrayList();
		
		// build the subsections
		int sectIndex = 0;
		for (FaultSection parentSect : fsd) {
			double ddw = parentSect.getOrigDownDipWidth();
			double maxSectLength = ddw*maxSubSectionLength;
			// the "2" here sets a minimum number of sub sections
			List<? extends FaultSection> newSubSects = parentSect.getSubSectionsList(maxSectLength, sectIndex, 2);
			subSections.addAll(newSubSects);
			sectIndex += newSubSects.size();
		}
		
		System.out.println(subSections.size()+" Sub Sections");
		
		// write subsection data to file
		File subSectDataFile = new File(outputDir, "sub_sections.xml");
		Document doc = XMLUtils.createDocumentWithRoot();
		FaultSystemIO.fsDataToXML(doc.getRootElement(), FaultModels.XML_ELEMENT_NAME, null, null, subSections);
		XMLUtils.writeDocumentToFile(subSectDataFile, doc);
		
		// instantiate our laugh test filter
		UCERF3PlausibilityConfig laughTest = UCERF3PlausibilityConfig.getDefault();
		// you will have to disable our coulomb filter as it uses a data file specific to our subsections
		CoulombRates coulombRates;
		if (coulombFilter) {
			coulombRates = remapCoulombRates(subSections, fm);
		} else {
			laughTest.setCoulombFilter(null);
			coulombRates = null;
		}
		laughTest.setCoulombRates(coulombRates);
		
		// calculate distances between each subsection
		Map<IDPairing, Double> subSectionDistances = DeformationModelFetcher.calculateDistances(maxDistance, subSections);
		Map<IDPairing, Double> reversed = Maps.newHashMap();
		// now add the reverse distance
		for (IDPairing pair : subSectionDistances.keySet()) {
			IDPairing reverse = pair.getReversed();
			reversed.put(reverse, subSectionDistances.get(pair));
		}
		subSectionDistances.putAll(reversed);
		Map<IDPairing, Double> subSectionAzimuths = DeformationModelFetcher.getSubSectionAzimuthMap(
				subSectionDistances.keySet(), subSections);
		
		// this separates the sub sections into clusters which are all within maxDist of each other and builds ruptures
		// fault model and deformation model here are needed by InversionFaultSystemRuptSet later, just to create a rup set
		// zip file
		SectionConnectionStrategy connectionStrategy = new UCERF3SectionConnectionStrategy(
				laughTest.getMaxJumpDist(), coulombRates);
		SectionClusterList clusters = new SectionClusterList(
				connectionStrategy, laughTest, subSections, subSectionDistances, subSectionAzimuths);
		
		List<List<Integer>> ruptures = Lists.newArrayList();
		for (SectionCluster cluster : clusters) {
			ruptures.addAll(cluster.getSectionIndicesForRuptures());
		}
		
		System.out.println("Created "+ruptures.size()+" ruptures");
		
		// write rupture/subsection associations to file
		// format: rupID	sectID1,sectID2,sectID3,...,sectIDN
		File rupFile = new File(outputDir, "ruptures.txt");
		FileWriter fw = new FileWriter(rupFile);
		Joiner j = Joiner.on(",");
		for (int i=0; i<ruptures.size(); i++) {
			fw.write(i+"\t"+j.join(ruptures.get(i))+"\n");
		}
		fw.close();
		
		// build actual rupture set for magnitudes and such
		LogicTreeBranch branch = LogicTreeBranch.fromValues(fm, DeformationModels.GEOLOGIC,
				ScalingRelationships.SHAW_2009_MOD, SlipAlongRuptureModels.TAPERED);
		InversionFaultSystemRupSet rupSet = new InversionFaultSystemRupSet(branch, clusters, subSections);
		
		File zipFile = new File(outputDir, "rupSet.zip");
		FaultSystemIO.writeRupSet(rupSet, zipFile);
	}
	
	public static CoulombRates remapCoulombRates(List<? extends FaultSection> subSections, FaultModels fm) throws IOException {
		// original coulomb rates
		CoulombRates origRates = CoulombRates.loadUCERF3CoulombRates(fm);
		
		// now load the original subsections
		List<? extends FaultSection> origSubSects = new DeformationModelFetcher(
				fm, DeformationModels.GEOLOGIC, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR,
				InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE).getSubSectionList();
		
		Map<IDPairing, CoulombRatesRecord> rates = Maps.newHashMap();
		
		Map<FaultSection, FaultSection> remappings = Maps.newHashMap();
		for (FaultSection sect : origSubSects)
			remappings.put(sect, getRemappedSubSect(sect, subSections));
		
		for (IDPairing pair : origRates.keySet()) {
			FaultSection subSect1 = origSubSects.get(pair.getID1());
			FaultSection subSect2 = origSubSects.get(pair.getID2());
			
			FaultSection mapped1 = remappings.get(subSect1);
			FaultSection mapped2 = remappings.get(subSect2);
			
			if (mapped1 != null && mapped2 != null) {
				CoulombRatesRecord origRec = origRates.get(pair);
				rates.put(new IDPairing(mapped1.getSectionId(), mapped2.getSectionId()), origRec);
			}
		}
		
		return new CoulombRates(null, rates);
	}
	
	private static FaultSection getRemappedSubSect(FaultSection origSubSect,
			List<? extends FaultSection> newSubSects) {
		for (FaultSection newSect : newSubSects) {
			if (newSect.getParentSectionId() != origSubSect.getParentSectionId())
				continue;
			if (newSect.getSectionName().equals(origSubSect.getSectionName()))
				return newSect;
		}
		return null;
	}

}
