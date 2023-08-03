package scratch.kevin.simulators.ruptures.multifaultSeparate;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.eq.MagUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.parsers.RSQSimFileWriter;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class MultifaultSeparator {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5413.instance();
		File catalogDir = catalog.getCatalogDir();
		
		File outputDir = new File(catalogDir.getParentFile(), catalogDir.getName()+"_multifault_separate");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minMag = 7d;
		double minFractForInclusion = 0.2;
		
		CSVFile<String> mappingCSV = new CSVFile<>(false);
		mappingCSV.addLine("Original Rupture ID", "Original Magnitude", "Modified Full Rupture ID", "Num Sub-Events",
				"Sub-Event ID 1", "Sub Event Magnitude 1", "...", "...", "Sub-Event ID N", "Sub Event Magnitude N");
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).skipYears(5000).load();
		System.out.println("Loaded "+events.size()+" events");
		
		FaultSystemSolution solution = null;
		boolean writeSol = false;
		
		File solFile = new File(outputDir, "solution_m"+(float)minMag+"_fract"+(float)minFractForInclusion+".zip");
		if (solFile.exists()) {
			solution = FaultSystemSolution.load(solFile);
		} else {
			solution = RSQSimUtils.buildFaultSystemSolution(
					catalog.getSubSects(), catalog.getElements(), events, 0d, minFractForInclusion);
			writeSol = true;
		}
		
		FaultSystemRupSet rupSet = solution.getRupSet();
		System.out.println(rupSet.getNumRuptures()+" ruptures");
		
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(catalog.getSubSects());
		
		RuptureConnectionSearch rsConnSearch = new RuptureConnectionSearch(rupSet, distAzCalc,
				10000d, RuptureConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		
		if (!rupSet.hasModule(ClusterRuptures.class)) {
			System.out.println("Building ClusterRuptures for RSQSim");
			rupSet.addModule(ClusterRuptures.instance(rupSet, rsConnSearch, false));
			writeSol = true;
		}
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		
		if (writeSol)
			solution.write(solFile);
		
		RSQSimSubSectionMapper mapper = new RSQSimSubSectionMapper(catalog.getSubSects(), catalog.getElements(), minFractForInclusion);
		NamedFaults named = catalog.getFaultModel().getNamedFaults();
		
		RSQSimStateTransitionFileReader transReader = catalog.getTransitions();
		RSQSimFileWriter writer = new RSQSimFileWriter(outputDir, "filtered", false, true, transReader.getVersion());
		
		List<SimulatorElement> allPatches = catalog.getElements();
		
		int debugID = 102;
		boolean exitAfterDebug = false;
		
		int numMultifault = 0;
		int numWritten = 0;
		double currentTime = 10000*SimulatorUtils.SECONDS_PER_YEAR;
		double timeBufferInternal = 60*60*24; // 1 day
		double timeBufferExternal = 60*60*24*7; // 7 days
		double maxTransWrittenTime = 0d;
		for (int i=0; i<events.size(); i++) {
			RSQSimEvent event = events.get(i);
			ClusterRupture cRup = cRups.get(i);
			
			if (cRup.getTotalNumClusters() > 1) {
				// might be a multifault rupture
				Map<String, List<FaultSubsectionCluster>> indvRups = new HashMap<>();
				
				for (FaultSubsectionCluster cluster : cRup.getClustersIterable()) {
					String name = named.getFaultName(cluster.parentSectionID);
					if (name == null)
						name = "Parent"+cluster.parentSectionID;
					List<FaultSubsectionCluster> clusters = indvRups.get(name);
					if (clusters == null) {
						clusters = new ArrayList<>();
						indvRups.put(name, clusters);
					}
					clusters.add(cluster);
				}
				
				if (indvRups.size() > 1) {
					System.out.println(event.getID()+" is a multifault rupture ("+indvRups.size()+" parts):\n\t"+cRup);
					numMultifault++;
					
					List<RSQSimStateTime> transitions = new ArrayList<>();
					transReader.getTransitions(event, transitions);
					
					int outputEventID = numMultifault*100;
					System.out.println("\tFull outputID="+outputEventID);
					System.out.println("\tOutput time="+(float)currentTime+" s = "+(float)(currentTime/SimulatorUtils.SECONDS_PER_YEAR)+" yrs");
					
					// first write full event
					double origTime = event.getTime();
					writer.writeTransitions(event, outputEventID, currentTime - origTime, transReader.getVersion(), transitions);
					RSQSimEvent eventModTime = event.cloneNewTime(currentTime, outputEventID);
					writer.writeEvent(eventModTime);
					
					List<String> line = new ArrayList<>(3 + 2*indvRups.size());
					line.add(event.getID()+"");
					line.add(event.getMagnitude()+"");
					line.add(outputEventID+"");
					line.add(indvRups.size()+"");
					
					if (outputEventID == debugID && exitAfterDebug)
						System.exit(1);
					
					outputEventID++;
					numWritten++;
					
					// now filtered versions
					for (String name : indvRups.keySet()) {
						List<FaultSubsectionCluster> clusters = indvRups.get(name);
						System.out.println("\tSub-Event "+outputEventID+": "+name+" => "+clusters);
						HashSet<Integer> includePatches = new HashSet<>();
						for (FaultSubsectionCluster cluster : clusters)
							for (FaultSection sect : cluster.subSects)
								for (SimulatorElement elem : mapper.getElementsForSection(sect))
									includePatches.add(elem.getID());
						
						// need to filter out the slip on other patches
						List<RSQSimEventRecord> modRecords = new ArrayList<>();
						for (EventRecord origRec : event) {
							Preconditions.checkState(origRec instanceof RSQSimEventRecord);
							RSQSimEventRecord rsRec = (RSQSimEventRecord)origRec;
							int[] recIDs = rsRec.getElementIDs();
							int numMatches = 0;
							for (int patch : recIDs)
								if (includePatches.contains(patch))
									numMatches++;
							
							if (numMatches == recIDs.length) {
								// copy directly
								modRecords.add(rsRec);
							} else if (numMatches > 0) {
								// filter it
								double[] recSlips = rsRec.getElementSlips();
								RSQSimEventRecord modRec = new RSQSimEventRecord(allPatches);
								
								double moment = 0d;
								
								for (int j=0; j<recIDs.length; j++) {
									int patchID = recIDs[j];
									if (includePatches.contains(patchID)) {
										double slip = recSlips[j];
										modRec.addSlip(patchID, slip);
										SimulatorElement elem = allPatches.get(patchID);
										Preconditions.checkState(elem.getID() == patchID);
										double patchArea = elem.getArea();
										moment += FaultMomentCalc.getMoment(patchArea, slip);
									}
								}
								modRec.setMoment(moment);
								modRecords.add(modRec);
							} // else skip
						}
						Preconditions.checkState(!modRecords.isEmpty());
						double sumMoment = 0d;
						for (RSQSimEventRecord rec : modRecords)
							sumMoment += rec.getMoment();
						double subMag = MagUtils.momentToMag(sumMoment);
						for (RSQSimEventRecord rec : modRecords)
							rec.setMagnitude(subMag);
						RSQSimEvent modEvent = new RSQSimEvent(modRecords);
						
						// filter transitions as well
						List<RSQSimStateTime> filteredTrans = new ArrayList<>();
						double firstRelTime = Double.POSITIVE_INFINITY;
						for (RSQSimStateTime trans : transitions) {
							if (includePatches.contains(trans.patchID)) {
								firstRelTime = Math.min(firstRelTime, trans.relativeTime);
								filteredTrans.add(trans);
							}
						}
						Preconditions.checkState(firstRelTime >= 0, "Bad firstRelTime=%s", firstRelTime);
						// now reset relative times
						currentTime += timeBufferInternal;
						System.out.println("\tSub-event time: "+currentTime);
						
						boolean debugTrans = debugID == outputEventID;
						for (int t=0; t<filteredTrans.size(); t++) {
							RSQSimStateTime trans = filteredTrans.get(t);
							double modRelTime = trans.relativeTime - firstRelTime;
							if (modRelTime < 0d) {
								Preconditions.checkState(modRelTime > -1e-16);
								modRelTime = 0d;
							}
							double modAbsTransTime = currentTime+modRelTime;
							Preconditions.checkState((float)modAbsTransTime >= (float)maxTransWrittenTime,
									"Transitions being written out of order? currentTime=%s, firstRelTime=%s, trans.relativeTime=%s, modRelTime=%s",
									currentTime, firstRelTime, trans.relativeTime, modRelTime);
							Preconditions.checkState(modAbsTransTime >= currentTime, "Trans before event time? %s vs %s",
									modAbsTransTime, modAbsTransTime);
							maxTransWrittenTime = modAbsTransTime;
							
							RSQSimStateTime modTrans = new RSQSimStateTime(modAbsTransTime,
									Float.max(0f, (float)modRelTime), outputEventID, trans.patchID, trans.state, trans.velocity);
							filteredTrans.set(t, modTrans);
							if (debugTrans)
								System.out.println("\tTransition "+t+"\n\t\tOriginal: "+trans+"\n\t\tModified: "+modTrans);
						}
						// time offset in transitions taken care of above
						writer.writeTransitions(modEvent, outputEventID, 0d, transReader.getVersion(), filteredTrans);
						modEvent = modEvent.cloneNewTime(currentTime, outputEventID);
						double modMoment = 0d;
						
						writer.writeEvent(modEvent, outputEventID);
						line.add(outputEventID+"");
						line.add(modEvent.getMagnitude()+"");
						
						if (outputEventID == debugID && exitAfterDebug)
							System.exit(1);
						
						outputEventID++;
						numWritten++;
					}
					mappingCSV.addLine(line);
					currentTime += timeBufferExternal;
				}
			}
		}
		
		writer.close();
		File paramFile = new File(catalogDir, "multiparam.in");
		File elemFile = new File(catalogDir, "zfault_Deepen.in");
		Files.copy(paramFile, new File(outputDir, paramFile.getName()));
		Files.copy(elemFile, new File(outputDir, elemFile.getName()));
		System.out.println("Found "+numMultifault+" multifault ruptures. Wrote "+numWritten+" parts");
		mappingCSV.writeToFile(new File(outputDir, "event_mappings.csv"));
	}

}
