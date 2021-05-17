package scratch.kevin.simulators.dists;

import java.util.Collections;
import java.util.List;

import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.EventsInWindowsMatcher;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.catBuild.CatalogBuilder;
import scratch.kevin.simulators.catBuild.ProbabalisticCatalogBuilder;
import scratch.kevin.simulators.catBuild.StandardCatalogBuilder;
import scratch.kevin.simulators.catBuild.StateBasedCatalogBuilder;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public enum RandomDistType {
		NORMAL("Random Normal Dist", "rand_norm_dist") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new NormalDistReturnPeriodProvider(rps);
			}
		},
		EXPONENTIAL("Random Exponential Dist", "rand_exp_dist") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new ExponentialDistReturnPeriodProvider(rps);
			}
		},
		LOG_NORMAL("Random Log-Normal Dist", "rand_lognorm_dist") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new LogNormalDistReturnPeriodProvider(rps);
			}
		},
		POISSON("Random Poisson Dist", "rand_poisson_dist") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new PoissonDistReturnPeriodProvider(rps);
			}
		},
		ACTUAL("Random Actual Dist", "rand_actual_dist") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new ActualDistReturnPeriodProvider(rps);
			}
		},
		STATE_BASED("N Dim. State Based", "rand_state_based") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				return new ActualDistReturnPeriodProvider(rps);
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new StateBasedCatalogBuilder();
			}
		},
		PREFERRED_SYN("Random Preferred Synthetic", "rand_preferred") {
			@Override
			public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				if (rupIden != null && rupIden instanceof ElementMagRangeDescription) {
					List<Integer> ids = ((ElementMagRangeDescription)rupIden).getElementIDs();
					if (ids.size() == 1 && ids.get(0) == ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID) {
						// San Jacinto
						List<RandomReturnPeriodProvider> provs = Lists.newArrayList();
						List<Double> weights = Lists.newArrayList();
						
						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(98), 0.2));
						weights.add(1.0);
						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(160), 0.15));
						weights.add(1.65);
						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(235), 0.18));
						weights.add(0.95);
						
//						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(98), 0.2));
//						weights.add(1.0);
//						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(160), 0.15));
//						weights.add(1.65);
//						provs.add(new LogNormalDistReturnPeriodProvider(Math.log(235), 0.1));
//						weights.add(1.0);
						
						CompoundDistReturnPeriodProvider comp = new CompoundDistReturnPeriodProvider(provs, weights);
						
						// now invert a bit
//						return RIDistPlot.invertForSJDist(RIDistPlot.getHistFunc(rps), comp, 1000, 0.3, 5);
						return comp;
					}
				}
				// default: log normal
				return LOG_NORMAL.instance(rupIden, rps, events);
			}
		},
		MOJAVE_DRIVER("Mojave Driver Dist", "rand_mojave_driver_dist") {
			private List<? extends SimulatorEvent> mojaveMatches;
			ElementMagRangeDescription mojaveIden;
			
			@Override
			public synchronized RandomReturnPeriodProvider instance(
					RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				Preconditions.checkState(rupIden instanceof ElementMagRangeDescription);
				ElementMagRangeDescription elemIden = (ElementMagRangeDescription)rupIden;
				if (mojaveMatches == null) {
					mojaveIden = new ElementMagRangeDescription(ElementMagRangeDescription.smartName("SAF Mojave", elemIden),
							ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, elemIden.getMinMag(), elemIden.getMaxMag());
					mojaveMatches = mojaveIden.getMatches(events);
				}
				if (elemIden.getElementIDs().contains(ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID))
					return ACTUAL.instance(rupIden, rps, events);
				
				List<? extends SimulatorEvent> followerMatches = rupIden.getMatches(events);
				return new FollowerReturnPeriodProvider(events, mojaveIden, mojaveMatches, rupIden, followerMatches, 10d, 1500);
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new ProbabalisticCatalogBuilder();
			}
		},
		MOJAVE_COACHELLA_CODRIVER("Mojave/Coachella Co-Driver Dist", "rand_coach_mojave_codriver_dist") {
			private List<? extends SimulatorEvent> mojaveMatches;
			ElementMagRangeDescription mojaveIden;
			private List<? extends SimulatorEvent> coachellaMatches;
			ElementMagRangeDescription coachellaIden;
			
			@Override
			public synchronized RandomReturnPeriodProvider instance(
					RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				Preconditions.checkState(rupIden instanceof ElementMagRangeDescription);
				ElementMagRangeDescription elemIden = (ElementMagRangeDescription)rupIden;
				if (mojaveMatches == null) {
					mojaveIden = new ElementMagRangeDescription(ElementMagRangeDescription.smartName("SAF Mojave", elemIden),
							ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, elemIden.getMinMag(), elemIden.getMaxMag());
					mojaveMatches = mojaveIden.getMatches(events);
				}
				if (coachellaMatches == null) {
					coachellaIden = new ElementMagRangeDescription(ElementMagRangeDescription.smartName("SAF Coachella", elemIden),
							ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, elemIden.getMinMag(), elemIden.getMaxMag());
					coachellaMatches = coachellaIden.getMatches(events);
				}
				if (elemIden.getElementIDs().contains(ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID))
//					return ACTUAL.instance(rupIden, rps, events);
					return new FollowerReturnPeriodProvider(events, coachellaIden, coachellaMatches, mojaveIden, mojaveMatches, 10d, 1500);
				
				List<? extends SimulatorEvent> followerMatches = rupIden.getMatches(events);
				return new FollowerReturnPeriodProvider(events, mojaveIden, mojaveMatches, rupIden, followerMatches, 10d, 1500);
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new ProbabalisticCatalogBuilder();
			}
		},
		SAN_JACINTO_COACHELLA_CODRIVER("San Jacinto/Coachella Co-Driver Dist", "rand_sj_coach_codriver_dist") {
			private List<? extends SimulatorEvent> sjMatches;
			ElementMagRangeDescription sjIden;
			private List<? extends SimulatorEvent> coachellaMatches;
			ElementMagRangeDescription coachellaIden;
			
			@Override
			public synchronized RandomReturnPeriodProvider instance(
					RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events) {
				Preconditions.checkState(rupIden instanceof ElementMagRangeDescription);
				ElementMagRangeDescription elemIden = (ElementMagRangeDescription)rupIden;
				if (sjMatches == null) {
					sjIden = new ElementMagRangeDescription(ElementMagRangeDescription.smartName("San Jacinto", elemIden),
							ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, elemIden.getMinMag(), elemIden.getMaxMag());
					sjMatches = sjIden.getMatches(events);
				}
				if (coachellaMatches == null) {
					coachellaIden = new ElementMagRangeDescription(ElementMagRangeDescription.smartName("SAF Coachella", elemIden),
							ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, elemIden.getMinMag(), elemIden.getMaxMag());
					coachellaMatches = coachellaIden.getMatches(events);
				}
				if (elemIden.getElementIDs().contains(ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID))
//					return ACTUAL.instance(rupIden, rps, events);
					return new FollowerReturnPeriodProvider(events, coachellaIden, coachellaMatches, sjIden, sjMatches, 10d, 1500);
				
				List<? extends SimulatorEvent> followerMatches = rupIden.getMatches(events);
				return new FollowerReturnPeriodProvider(events, sjIden, sjMatches, rupIden, followerMatches, 10d, 1500);
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new ProbabalisticCatalogBuilder();
			}
		},
		PROBABILISTIC_SHUFFLE("Prob. Shuffle", "rand_prob_shuffle") {
			@Override
			public RandomReturnPeriodProvider instance(
					RuptureIdentifier rupIden, double[] rps,
					List<? extends SimulatorEvent> events) {
				final double rate = (rps.length+1d)/(events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears());
				return new ProbabilisticReturnPeriodProvider() {
					
					@Override
					public double getReturnPeriod() {
						return Double.NaN;
					}
					
					@Override
					public double getPreferredWindowLength() {
						return 1d;
					}
					
					@Override
					public PossibleRupture getPossibleRupture(List<SimulatorEvent> prevEvents,
							double windowStart, double windowEnd) {
						double duration = windowEnd - windowStart;
						double prob = 1-Math.exp(-rate*duration);
						return new PossibleRupture(prob, windowStart+0.5*duration);
					}
				};
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new ProbabalisticCatalogBuilder();
			}
		},
		TRUE_SHUFFLE("True Shuffle", "rand_true_shuffle") {
			@Override
			public RandomReturnPeriodProvider instance(
					RuptureIdentifier rupIden, double[] rps,
					List<? extends SimulatorEvent> events) {
				return null;
			}

			@Override
			public CatalogBuilder getBuilder() {
				return new CatalogBuilder() {
					
					@Override
					public List<SimulatorEvent> buildCatalog(List<? extends SimulatorEvent> events,
							List<RandomReturnPeriodProvider> randomRPsList,
							List<List<? extends SimulatorEvent>> eventListsToResample, boolean trim) {
						
						double startTime = events.get(0).getTime();
						double endTime = events.get(events.size()-1).getTime();
						
						double duration = endTime - startTime;
						
						List<SimulatorEvent> randEvents = Lists.newArrayList();
						for (SimulatorEvent event : events) {
							double newTime = Math.random()*duration + startTime;
							randEvents.add(event.cloneNewTime( newTime, event.getID()));
						}
						Collections.sort(randEvents);
						int eventID = 0;
						for (SimulatorEvent e : randEvents)
							e.setID(eventID++);
						
						return randEvents;
					}
				};
			}
		};
		
		private String name, fNameAdd;
		private RandomDistType(String name, String fNameAdd) {
			this.name = name;
			this.fNameAdd = fNameAdd;
		}
		public String getName() {
			return name;
		}
		public String getFNameAdd() {
			return fNameAdd;
		}
		public abstract RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, double[] rps, List<? extends SimulatorEvent> events);
		
		public RandomReturnPeriodProvider instance(RuptureIdentifier rupIden, List<? extends SimulatorEvent> events) {
			List<? extends SimulatorEvent> matches = rupIden.getMatches(events);
			double[] rps = PeriodicityPlotter.getRPs(matches);
			return instance(rupIden, rps, events);
		}
		
		public CatalogBuilder getBuilder() {
			return new StandardCatalogBuilder();
		}
	}