package scratch.kevin.simulators.catBuild;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.EventsInWindowsMatcher;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.kevin.simulators.dists.FollowerReturnPeriodProvider;
import scratch.kevin.simulators.dists.PossibleRupture;
import scratch.kevin.simulators.dists.ProbabilisticReturnPeriodProvider;
import scratch.kevin.simulators.dists.RandomReturnPeriodProvider;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

public class ProbabalisticCatalogBuilder implements CatalogBuilder {
		
		private StandardCatalogBuilder standardBuild = new StandardCatalogBuilder();
		private Random rand = new Random();
		
		private static final boolean D = true;

		@Override
		public List<SimulatorEvent> buildCatalog(List<? extends SimulatorEvent> events,
				List<RandomReturnPeriodProvider> randomRPsList,
				List<List<? extends SimulatorEvent>> eventListsToResample, boolean trim) {
			
			// separate standard from regular
			List<RandomReturnPeriodProvider> standardRPs = Lists.newArrayList();
			List<List<? extends SimulatorEvent>> standardEventLists = Lists.newArrayList();
			List<ProbabilisticReturnPeriodProvider> probRPs = Lists.newArrayList();
			List<List<? extends SimulatorEvent>> probEventLists = Lists.newArrayList();
			for (int i=0; i<randomRPsList.size(); i++) {
				RandomReturnPeriodProvider rp = randomRPsList.get(i);
				if (rp instanceof ProbabilisticReturnPeriodProvider) {
					probRPs.add((ProbabilisticReturnPeriodProvider)rp);
					List<SimulatorEvent> eventList = Lists.newArrayList(eventListsToResample.get(i));
					Collections.shuffle(eventList);
					probEventLists.add(eventList);
				} else {
					standardRPs.add(rp);
					standardEventLists.add(eventListsToResample.get(i));
				}
			}
			
			// first populate any standard ones
			List<SimulatorEvent> standardEvents;
			if (!standardRPs.isEmpty())
				standardEvents = standardBuild.buildCatalog(events, standardRPs, standardEventLists, true);
			else
				standardEvents = Lists.newArrayList();
			
			if (probRPs.isEmpty())
				return standardEvents;
			
			double startTime = 0;
			double maxTime = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
			double timeDelta = probRPs.get(0).getPreferredWindowLength();
			// make sure all time deltas are the same (for now at least)
			// TODO allow variable deltas
			for (ProbabilisticReturnPeriodProvider probRP : probRPs)
				Preconditions.checkState((float)timeDelta == (float)probRP.getPreferredWindowLength());
			double halfDelta = timeDelta*0.5;
			
			int standardEventIndex = 0;
			List<SimulatorEvent> runningEvents = Lists.newArrayList();
			int probAdded = 0;
			
			int[] probEventIndexes = new int[probRPs.size()];
			
			// do a for loop to avoid propagating double precision errors
			int numSteps = (int)((maxTime - startTime)/timeDelta);
			System.out.println("Doing probabilisic build with "+numSteps+" steps and "+probRPs.size()+" providers!");
			System.out.println("Already have "+standardEvents.size()+" standard events");
			
			Stopwatch loopWatch = null;
			Stopwatch probWatch = null;
			if (D) {
				loopWatch =Stopwatch.createStarted();
				probWatch = Stopwatch.createUnstarted();
			}
			
			List<Integer> iterationIndexes = Lists.newArrayList();
			for (int i=0; i<probRPs.size(); i++)
				iterationIndexes.add(i);
			
//			int stepDiscr = 100;
//			double stepProbMult = 1d/(double)stepDiscr;
//			
//			double subStepDelta = timeDelta/(double)stepDiscr*0.5;
//			EvenlyDiscretizedFunc stepDiscrTimes = new EvenlyDiscretizedFunc(subStepDelta, timeDelta-subStepDelta, stepDiscr);
			
			int eventID = 0;
			
			for (int step=0; step<numSteps; step++) {
				double windowStart = startTime + timeDelta*step;
				double windowEnd = windowStart + timeDelta;
				
				if (D && step % 10000 == 0) {
					double time = startTime + halfDelta + timeDelta*step;
					long loopSecs = loopWatch.elapsed(TimeUnit.SECONDS);
					long probSecs = probWatch.elapsed(TimeUnit.SECONDS);
					double probFract = (double)probSecs/(double)loopSecs;
					double timeStepMillis = (double)loopWatch.elapsed(TimeUnit.MILLISECONDS)/(double)step;
					System.out.println("Step "+step+",\ttime="+(float)time+"\tevents="+runningEvents.size()
							+"\tprobEvents="+probAdded+"\tloopSecs="+loopSecs+",\tprobSecs="+probSecs
							+"\tprobFract="+(float)probFract+"\tstepMillis="+(float)timeStepMillis);
				}
					
				// populate any standard events before the current time
				for (int i=standardEventIndex; i<standardEvents.size(); i++) {
					SimulatorEvent e = standardEvents.get(i);
					if (e.getTimeInYears()<=windowEnd) {
//						runningEvents.add(e);
						SimulatorEvent newE = e.cloneNewTime(0.5*(windowStart+windowEnd)*General_EQSIM_Tools.SECONDS_PER_YEAR, eventID++);
						runningEvents.add(newE);
						standardEventIndex = i+1;
					} else {
						break;
					}
				}
				
				List<SimulatorEvent> eventsToAdd = Lists.newArrayList();

				// now do probabilistic ones
				
				Collections.shuffle(iterationIndexes);
				for (int i:iterationIndexes) {
					ProbabilisticReturnPeriodProvider probRP = probRPs.get(i);

					if (D) probWatch.start();
					PossibleRupture rup = probRP.getPossibleRupture(runningEvents, windowStart, windowEnd);
					if (D) probWatch.stop();
					if (rup == null)
						continue;
					double prob = rup.getProb();
					if (Double.isNaN(prob))
						continue;
					double r = rand.nextDouble();
					if (r<prob) {
						// add an event
						List<? extends SimulatorEvent> myEvents = probEventLists.get(i);

						int ind = probEventIndexes[i]++;
						if (ind == myEvents.size()) {
							// roll back to start
							ind = 0;
							probEventIndexes[i] = 0;
						}

						SimulatorEvent e = myEvents.get(ind);
						double rupTime = rup.getEventTimeYears();
//						if (rupTime > windowEnd) {
//							// see if any standard events in this little delta
//							// TODO remove debug
//							for (int j=standardEventIndex; j<standardEvents.size(); j++) {
//								double newTime = standardEvents.get(j).getTimeInYears();
//								if (newTime<=rupTime) {
//									System.out.println("Would insert an extra before rup time!!! window: "
//											+(float)windowStart+"=>"+(float)windowEnd+", rupTime="+(float)rupTime
//											+", newTime="+(float)newTime);
//								} else {
//									break;
//								}
//							}
//						}
						double timeSecs = rupTime * General_EQSIM_Tools.SECONDS_PER_YEAR;
						SimulatorEvent newE = e.cloneNewTime(timeSecs, eventID++);

						boolean inserted = false;
						for (int checkIndex=runningEvents.size(); --checkIndex>=0;) {
							if (runningEvents.get(checkIndex).getTime() <= timeSecs) {
								runningEvents.add(checkIndex+1, newE);
								inserted = true;
								break;
							}
						}
						if (!inserted)
							runningEvents.add(0, newE);
//						runningEvents.add(newE);
						probAdded++;
					}
				}
				
//				Collections.sort(eventsToAdd);
//				runningEvents.addAll(eventsToAdd);
			}
			
			for (int i=0; i<probRPs.size(); i++) {
				if (probRPs.get(i) instanceof FollowerReturnPeriodProvider) {
					System.out.println("\nFollower Stats:");
					((FollowerReturnPeriodProvider)probRPs.get(i)).printStats();
				}
			}
			
			Collections.sort(runningEvents);
			
			return runningEvents;
		}
		
	}