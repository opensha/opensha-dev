package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SectionIDIden;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import java.util.List;
import java.util.Map;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimMarkovChainBuilder {
	
	public static List<RuptureIdentifier> getParentFaultIdens(RSQSimCatalog catalog, boolean middleSubSect, String... parentNames) {
		return getParentFaultIdens(catalog, middleSubSect, Lists.newArrayList(parentNames));
	}
	
	public static List<RuptureIdentifier> getParentFaultIdens(RSQSimCatalog catalog, boolean middleSubSect, List<String> parentNames) {
		List<RuptureIdentifier> faultIdens = new ArrayList<>();
		
		List<FaultSectionPrefData> u3SubSects = catalog.getU3SubSects();
		Map<String, List<FaultSectionPrefData>> parentNameToSubSectsMap = new HashMap<>();
		for (FaultSectionPrefData sect : u3SubSects) {
			String parent = sect.getParentSectionName();
			List<FaultSectionPrefData> parentSects = parentNameToSubSectsMap.get(parent);
			if (parentSects == null) {
				parentSects = new ArrayList<>();
				parentNameToSubSectsMap.put(parent, parentSects);
			}
			parentSects.add(sect);
		}
		
		List<SimulatorElement> elems;
		try {
			elems = catalog.getElements();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		for (String parent : parentNames) {
			Preconditions.checkState(parentNameToSubSectsMap.containsKey(parent));
			List<FaultSectionPrefData> sects = parentNameToSubSectsMap.get(parent);
			
			List<Integer> sectIDs = new ArrayList<>();
			for (FaultSectionPrefData sect : sects)
				sectIDs.add(sect.getSectionId());
			Collections.sort(sectIDs);
			
			SectionIDIden iden;
			if (middleSubSect)
				iden = new SectionIDIden(parent, elems, sectIDs.get(sectIDs.size()/2));
			else
				iden = new SectionIDIden(parent, elems, sectIDs);
			
			faultIdens.add(iden);
		}
		
		return faultIdens;
	}
	
	public static void main(String[] args) throws IOException {
		File catalogsDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(catalogsDir);
		
		List<String> parentNames = new ArrayList<>();
		
		parentNames.add("San Andreas (Carrizo) rev");
		parentNames.add("San Andreas (Cholame) rev");
		parentNames.add("San Andreas (Mojave S)");
		parentNames.add("San Andreas (Coachella) rev");
		parentNames.add("San Jacinto (Anza) rev");
		parentNames.add("Garlock (West)");
		
		double minMag = 7d;
		double distSpacing = 10d; // years
		double skipYears = 5000;
		boolean middleSubSect = true; // else any
		
		List<RuptureIdentifier> faultIdens = getParentFaultIdens(catalog, middleSubSect, parentNames);
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).matches(new LogicalOrRupIden(faultIdens)).load();
		
		List<int[]> path = MarkovChainBuilder.getStatesPath(distSpacing, events, faultIdens, skipYears);
		System.out.println("Writing Markov Chain with "+path.size()+" states");
		
		File outputFile = new File(catalog.getCatalogDir(), "markov_chain_m"+(float)minMag+"_"+faultIdens.size()+"sects.csv");
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Index");
		for (RuptureIdentifier iden : faultIdens)
			header.add(iden.getName());
		csv.addLine(header);
		
		for (int i=0; i<path.size(); i++) {
			List<String> line = new ArrayList<>();
			line.add(i+"");
			for (int state : path.get(i))
				line.add(state+"");
			csv.addLine(line);
		}
		
		System.out.println("Writing to: "+outputFile.getAbsolutePath());
		csv.writeToFile(outputFile);
	}

}
