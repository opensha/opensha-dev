package org.opensha.nshmp.shaded.model;

import static java.util.stream.Collectors.toList;
import static org.opensha.sha.util.TectonicRegionType.ACTIVE_SHALLOW;
import static org.opensha.sha.util.TectonicRegionType.STABLE_SHALLOW;
import static org.opensha.sha.util.TectonicRegionType.SUBDUCTION_INTERFACE;
import static org.opensha.sha.util.TectonicRegionType.SUBDUCTION_SLAB;
import static org.opensha.sha.util.TectonicRegionType.VOLCANIC;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.numbers.core.Precision;
import org.opensha.commons.data.TimeSpan;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

import org.opensha.nshmp.shaded.data.NshmpIndexing;
import org.opensha.nshmp.shaded.gmm.NshmpGmm;
import org.opensha.nshmp.shaded.model.NshmpSystemRuptureSet.SystemRupture;
import org.opensha.nshmp.shaded.tree.NshmpBranch;

public class NshmErf extends AbstractERF {

	private final NshmpHazardModel model;
	private final List<NshmSource> allSources;
	private final Multimap<TectonicRegionType, NshmSource> sourceMap;

	private final boolean activeCrust;
	private final boolean stableCrust;
	private final boolean subInterface;
	private final boolean subSlab;
	private final boolean volcanic;
	private final boolean grid;
	private final boolean faults;
	private Set<TectonicRegionType> trts;

	public NshmErf(Path path, Set<TectonicRegionType> trts,
			IncludeBackgroundOption gridOption) {
		this(NshmpHazardModel.load(path), trts, gridOption);
	}

	public NshmErf(NshmpHazardModel model, Set<TectonicRegionType> trts,
			IncludeBackgroundOption gridOption) {
		this.model = model;
		this.trts = trts;
		allSources = new ArrayList<>();
		sourceMap = MultimapBuilder
				.enumKeys(TectonicRegionType.class)
				.arrayListValues()
				.build();

		activeCrust = trts.contains(ACTIVE_SHALLOW) || trts.isEmpty();
		stableCrust = trts.contains(STABLE_SHALLOW) || trts.isEmpty();
		subInterface = trts.contains(SUBDUCTION_INTERFACE) || trts.isEmpty();
		subSlab = trts.contains(SUBDUCTION_SLAB) || trts.isEmpty();
		volcanic = trts.contains(VOLCANIC) || trts.isEmpty();
		this.grid = gridOption == IncludeBackgroundOption.INCLUDE ||
				gridOption == IncludeBackgroundOption.ONLY;
		this.faults = gridOption != IncludeBackgroundOption.ONLY;

		init();
	}
	
	private static final boolean PROCESS_TRT_OVERRIDES = true;

	private void init() {

		// ERF initializers
		timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		timeSpan.setDuration(50.0);
		setTimeSpan(timeSpan);

		// nshmp-haz initializers
		Multimap<NshmpTectonicSetting, NshmpSourceTree> trees = model.trees();
		if (PROCESS_TRT_OVERRIDES) {
			List<NshmpTRTBranch> trtMappedRuptureSets = parseTRTOverrides(trees);
			for (NshmpTRTBranch trtBranch : trtMappedRuptureSets) {
				if (!trts.isEmpty() && !trts.contains(trtBranch.trt))
					continue;
				List<NshmSource> sources = initBranch(trtBranch.branch);
				sources.forEach(s -> s.setTectonicRegionType(trtBranch.trt));
				allSources.addAll(sources);
				sourceMap.putAll(trtBranch.trt, sources);
			}
			allSources.sort(new Comparator<NshmSource>() {
				@Override
				public int compare(NshmSource o1, NshmSource o2) {
					return Integer.compare(o1.getNSHM_ID(), o2.getNSHM_ID());
				}
			});
		} else {
			for (Entry<NshmpTectonicSetting, NshmpSourceTree> entry : trees.entries()) {
				NshmpTectonicSetting setting = entry.getKey();
				NshmpSourceTree tree = entry.getValue();
				NshmpSourceType type = tree.type();

				TectonicRegionType trt = NshmUtil.tectonicSettingToType(setting, type);
				if (!trts.isEmpty() && !trts.contains(trt))
					continue;
				List<NshmSource> sources = initTree(tree);
				sources.forEach(s -> s.setTectonicRegionType(trt));
				allSources.addAll(sources);
				sourceMap.putAll(trt, sources);
			}
		}
	}
	
	private record NshmpTRTBranch(TectonicRegionType trt, NshmpBranch<NshmpRuptureSet> branch) {}
	
	/**
	 * Map each branch to a tectonic region type, also processing any GMM-overrides that Peter might sneak in there
	 * to change the TRT within a different NshmpTectonicSetting
	 * @param trees
	 * @return
	 */
	private static List<NshmpTRTBranch> parseTRTOverrides(
			Multimap<NshmpTectonicSetting, NshmpSourceTree> trees) {
		List<NshmpTRTBranch> ret = new ArrayList<>();
		
		for (Entry<NshmpTectonicSetting, NshmpSourceTree> entry : trees.entries()) {
			NshmpTectonicSetting setting = entry.getKey();
			NshmpSourceTree tree = entry.getValue();
			NshmpSourceType type = tree.type();
			TectonicRegionType origTRT = NshmUtil.tectonicSettingToType(setting, type);
			
			System.out.println("Processing tree "+tree.id()+". "+tree.name()+" ("+origTRT.name()+")");
			
			for (NshmpBranch<NshmpRuptureSet> branch : tree) {
				if (origTRT == VOLCANIC) {
					// copy over as is, no current GMM type for volcanic
					ret.add(new NshmpTRTBranch(origTRT, branch));
					continue;
				}
				NshmpRuptureSet rs = branch.value();
				System.out.println("\tProcessing RS "+rs.id()+". "+rs.name()+" ("+origTRT.name()+")");
				NshmpGmmTree gmmTree = rs.gmmTree();
				double weightOrig = 0d;
				double weightOverrideSum = 0d;
				Map<TectonicRegionType, Double> overrideWeights = new EnumMap<>(TectonicRegionType.class);
				for (NshmpBranch<NshmpGmm> gmm : gmmTree.tree()) {
					TectonicRegionType gmmTRT = NSHMP_GMM_Wrapper.trtForType(gmm.value().type());
					if (gmmTRT != origTRT) {
						// we have a GMM TRT override
						weightOverrideSum += gmm.weight();
						if (overrideWeights.containsKey(gmmTRT))
							overrideWeights.put(gmmTRT, overrideWeights.get(gmmTRT)+gmm.weight());
						else
							overrideWeights.put(gmmTRT, gmm.weight());
					} else {
						weightOrig += gmm.weight();
					}
				}
				if (weightOverrideSum > 0d) {
					// we have overrides to process
					Preconditions.checkState(Precision.equals(weightOverrideSum + weightOrig, 1d, 1e-4));
					System.err.println("Detected a GMM TRT override for source "+rs.id()+". "+rs.name()+" with original TRT="+origTRT.name());
					if (weightOrig > 0d) {
						System.err.println("\t"+origTRT.name()+":\t"+weightOrig);
						ret.add(new NshmpTRTBranch(origTRT, new WeightScaledBranch(branch, weightOrig, origTRT)));
					}
					for (TectonicRegionType trt : overrideWeights.keySet()) {
						double weight = overrideWeights.get(trt);
						System.err.println("\t"+trt.name()+":\t"+weight);
						ret.add(new NshmpTRTBranch(trt, new WeightScaledBranch(branch, weight, trt)));
					}
				} else {
					// copy over as is
					ret.add(new NshmpTRTBranch(origTRT, branch));
				}
			}
		}
		
		return ret;
	}
	
	private static class WeightScaledBranch implements NshmpBranch<NshmpRuptureSet> {
		
		private NshmpBranch<NshmpRuptureSet> upstream;
		private double scale;
		private TectonicRegionType trt;

		public WeightScaledBranch(NshmpBranch<NshmpRuptureSet> upstream, double scale, TectonicRegionType trt) {
			this.upstream = upstream;
			this.scale = scale;
			this.trt = trt;
		}

		@Override
		public String id() {
			return upstream.id()+"-"+trt.name();
		}

		@Override
		public NshmpRuptureSet value() {
			return upstream.value();
		}

		@Override
		public double weight() {
			return upstream.weight()*scale;
		}
		
	}

	public List<NshmSource> allSources() {
		return allSources;
	}

	public Multimap<TectonicRegionType, NshmSource> sourceMap() {
		return sourceMap;
	}

	private List<NshmSource> initTree(NshmpSourceTree tree) {
		List<NshmSource> sources = new ArrayList<>();
		double duration = getTimeSpan().getDuration();
		tree.stream()
		.map(branch -> sourcesFromBranch(branch, duration))
		.forEach(sources::addAll);

		// tree.stream()
		// .map(branch -> {
		// List<NshmSource> brSrcs = sourcesFromBranch(branch, duration);
		// if (brSrcs.size() > 0) {
		// System.out.println("type: " + branch.value().type());
		// }
		// return brSrcs;
		// })
		// .forEach(list -> {
		//// if (list.size() > 0) {
		//// System.out.println("br: " + list.get(0).getTectonicRegionType());
		//// }
		// sources.addAll(list);
		// });

		sources.sort(new Comparator<NshmSource>() {
			@Override
			public int compare(NshmSource o1, NshmSource o2) {
				return Integer.compare(o1.getNSHM_ID(), o2.getNSHM_ID());
			}
		});

		return sources;
	}

	private List<NshmSource> initBranch(NshmpBranch<NshmpRuptureSet> branch) {
		List<NshmSource> sources = new ArrayList<>();
		double duration = getTimeSpan().getDuration();
		sources.addAll(sourcesFromBranch(branch, duration));

		return sources;
	}

	private List<NshmSource> sourcesFromBranch(
			NshmpBranch<NshmpRuptureSet> branch,
			double duration) {

		NshmpRuptureSet ruptureSet = branch.value();
		double weight = branch.weight();

		switch (ruptureSet.type()) {

		case GRID:
			NshmpGriddedRuptureSet grs = (NshmpGriddedRuptureSet) ruptureSet;
			return (grid)
					? gridRuptureSetToSources(grs, weight, duration)
							: List.of();

		case ZONE:
			NshmpGriddedRuptureSet zrs = (NshmpGriddedRuptureSet) ruptureSet;
			return (faults)
					? gridRuptureSetToSources(zrs, weight, duration)
							: List.of();

		case INTERFACE_GRID:
			NshmpGriddedRuptureSet igrs = (NshmpGriddedRuptureSet) ruptureSet;
			return (subInterface && grid)
					? gridRuptureSetToSources(igrs, weight, duration)
							: List.of();

		case INTRASLAB_GRID:
			NshmpGriddedRuptureSet isgrs = (NshmpGriddedRuptureSet) ruptureSet;
			return (subSlab && grid)
					? gridRuptureSetToSources(isgrs, weight, duration)
							: List.of();

		case SLAB:
			NshmpGriddedRuptureSet slabRuptures = (NshmpGriddedRuptureSet) ruptureSet;
			return (subSlab && grid)
					? gridRuptureSetToSources(slabRuptures, weight, duration)
							: List.of();

		case FAULT_CLUSTER:
			NshmpClusterRuptureSet crs = (NshmpClusterRuptureSet) ruptureSet;
			return (faults)
					? clusterRuptureSetToSources(crs, weight, duration)
							: List.of();

		case FAULT_SYSTEM:
			NshmpSystemRuptureSet srs = (NshmpSystemRuptureSet) ruptureSet;
			return (faults)
					? systemRuptureSetToSources(srs, weight, duration)
							: List.of();

		case INTERFACE:
			return (subInterface && faults)
					? iterableRuptureSetToSources((NshmpIterableRuptureSet) ruptureSet, weight, duration)
							: List.of();

		case INTERFACE_CLUSTER:
			NshmpClusterRuptureSet icrs = (NshmpClusterRuptureSet) ruptureSet;
			return (subInterface && faults)
					? clusterRuptureSetToSources(icrs, weight, duration)
							: List.of();

		case INTERFACE_SYSTEM:
			NshmpSystemRuptureSet isrs = (NshmpSystemRuptureSet) ruptureSet;
			return (faults)
					? systemRuptureSetToSources(isrs, weight, duration)
							: List.of();

		default:
			return (faults)
					? iterableRuptureSetToSources((NshmpIterableRuptureSet) ruptureSet, weight, duration)
							: List.of();
		}
	}

	private static List<NshmSource> gridRuptureSetToSources(
			NshmpGriddedRuptureSet ruptureSet,
			double weight,
			double duration) {

		List<NshmSource> sources = new ArrayList<>();
		for (NshmpGridSource gridSource : ruptureSet) {
			sources.add(new NshmSource.Point(gridSource, weight, duration));
		}
		return sources;
	}

	private static List<NshmSource> iterableRuptureSetToSources(
			NshmpIterableRuptureSet ruptureSet,
			double weight,
			double duration) {

		return List.of(new NshmSource.Fault(ruptureSet, weight, duration));
	}

	private static List<NshmSource> systemRuptureSetToSources(
			NshmpSystemRuptureSet srs,
			double weight,
			double duration) {

		List<NshmSurface> surfaces = Arrays.stream(srs.sections)
				.map(section -> new NshmSurface(section))
				.collect(Collectors.toList());

		// NshmpSystemRuptureSet.stream() not supported but should be.
		// Iterator should work, even if it isn't used in nshm calc
		// pathways

		List<NshmSource> sources = new ArrayList<>(srs.size());
		for (int i = 0; i < srs.size(); i++) {
			SystemRupture source = (SystemRupture) srs.get(i);
			int[] sectionIndices = NshmpIndexing.bitsToIndices(source.bitset());
			List<NshmSurface> ruptureSurfaces = IntStream.of(sectionIndices)
					.mapToObj(surfaces::get)
					.collect(Collectors.toList());
			sources.add(new NshmSource.System(srs, source, weight, duration, ruptureSurfaces));
		}
		return sources;
	}

	private static List<NshmSource> clusterRuptureSetToSources(
			NshmpClusterRuptureSet crs,
			double weight,
			double duration) {

		double rate = crs.rate();
		return crs.ruptureSets().stream()
				.map(rs -> new NshmSource.Fault(rs, weight * rate, duration))
				.collect(toList());
	}

	@Override
	public int getNumSources() {
		return allSources.size();
	}

	@Override
	public ProbEqkSource getSource(int index) {
		return allSources.get(index);
	}

	@Override
	public void updateForecast() {
		double duration = getTimeSpan().getDuration();
		allSources.forEach(src -> src.setDuration(duration));
	}

	@Override
	public String getName() {
		return model.name();
	}

}
