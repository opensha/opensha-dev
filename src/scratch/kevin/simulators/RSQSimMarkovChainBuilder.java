package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SectionIDIden;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;

import java.util.List;
import java.util.Map;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimMarkovChainBuilder {
	
	public static SectionIDIden getU3_SectionIdentifier(RSQSimCatalog catalog, double minAreaFract, boolean middleSubSectOnly,
			String... u3ParentSectionNames) throws IOException {
		int[] parentIDs = new int[u3ParentSectionNames.length];
		for (int i=0; i<parentIDs.length; i++) {
			String name = u3ParentSectionNames[i];
			parentIDs[i] = -1;
			for (FaultSectionPrefData sect : catalog.getU3SubSects()) {
				if (sect.getParentSectionName().equals(name)) {
					parentIDs[i] = sect.getParentSectionId();
					break;
				}
			}
			Preconditions.checkState(parentIDs[i] >= 0, "Index not found for parent %s", name);
		}
		return getU3_SectionIdentifier(catalog, minAreaFract, middleSubSectOnly, parentIDs);
	}
	
	public static SectionIDIden getU3_SectionIdentifier(RSQSimCatalog catalog, double minAreaFract, boolean middleSubSectOnly,
			int... u3ParentSectionIDs) throws IOException {
		List<FaultSectionPrefData> subSects = catalog.getU3SubSects();
		// sometimes the RSQSim section numbers are 1-based inssead of 0-based, this is used to convert to UCERF3 convention
		int subSectOffset = RSQSimUtils.getSubSectIndexOffset(catalog.getElements(), subSects);
		
		String[] parentNames = new String[u3ParentSectionIDs.length];
		for (int i=0; i<u3ParentSectionIDs.length; i++) {
			for (FaultSectionPrefData sect : subSects) {
				if (sect.getParentSectionId() == u3ParentSectionIDs[i]) {
					parentNames[i] = sect.getParentSectionName();
					break;
				}
			}
			Preconditions.checkNotNull(parentNames[i]);
		}
		
		String parentName = getCombinedParentsName(subSects, parentNames);
		
		List<Integer> sectIDsForParent = new ArrayList<Integer>();
		for (FaultSectionPrefData subSect : subSects) {
			if (Ints.contains(u3ParentSectionIDs, subSect.getParentSectionId())) {
				// this subsection is on the given parent section
				sectIDsForParent.add(subSect.getSectionId()+subSectOffset);
			}
		}
		Preconditions.checkState(!sectIDsForParent.isEmpty(),
				"No subsections found with parent section ID="+Joiner.on(",").join(Ints.asList(u3ParentSectionIDs)));
		SectionIDIden iden;
		if (middleSubSectOnly)
			iden = new SectionIDIden(parentName, catalog.getElements(), sectIDsForParent.get(sectIDsForParent.size()/2));
		else
			iden = new SectionIDIden(parentName, catalog.getElements(), sectIDsForParent);
		iden.setAreaFractForInclusion(minAreaFract);
		return iden;
	}
	
	public static String getCombinedParentsName(List<FaultSectionPrefData> subSects, String[] parentNames) {
		String parentName;
		if (parentNames.length > 1) {
			String prefix = StringUtils.getCommonPrefix(parentNames);
			if (prefix.contains("("))
				prefix = prefix.substring(0, prefix.indexOf("("));
			List<String> suffixes = new ArrayList<>();
			if (prefix.length() > 5) {
				for (String name : parentNames) {
					String suffix = name.substring(prefix.length()).trim();
					if (suffix.startsWith("("))
						suffix = suffix.substring(1).trim();
					if (prefix.contains("(") && suffix.contains(")"));
						suffix = suffix.substring(0, suffix.indexOf(")"));
					suffixes.add(suffix);
				}
				parentName = prefix.trim()+" ["+Joiner.on("; ").join(suffixes)+"]";
			} else {
				parentName = Joiner.on("; ").join(Lists.newArrayList(parentNames));
			}
		} else {
			parentName = parentNames[0];
		}
		return parentName;
	}
	
	public static List<RuptureIdentifier> getU3_SectionIdentifiers(RSQSimCatalog catalog, double minAreaFract,
			boolean middleSubSect, String... parentNames) throws IOException {
		return getU3_SectionIdentifiers(catalog, minAreaFract, middleSubSect, Lists.newArrayList(parentNames));
	}
	
	public static List<RuptureIdentifier> getU3_SectionIdentifiers(RSQSimCatalog catalog, double minAreaFract,
			boolean middleSubSect, List<String> parentNames) throws IOException {
		List<RuptureIdentifier> faultIdens = new ArrayList<>();
		
		for (String parentName : parentNames)
			faultIdens.add(getU3_SectionIdentifier(catalog, minAreaFract, middleSubSect, parentName));
		
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
		boolean middleSubSect = false; // else any
		double minAreaFract = 0.2;
		
		List<RuptureIdentifier> faultIdens = getU3_SectionIdentifiers(catalog, minAreaFract, middleSubSect, parentNames);
		
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
