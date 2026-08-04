package gov.usgs.earthquake.nshmp.model;

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
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.opensha.commons.data.TimeSpan;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

import gov.usgs.earthquake.nshmp.data.Indexing;
import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemRupture;
import gov.usgs.earthquake.nshmp.tree.Branch;

public class NshmErf extends AbstractERF {

	private final HazardModel model;
	private final List<NshmSource> allSources;
	private final Multimap<TectonicRegionType, NshmSource> sourceMap;

	private final boolean activeCrust;
	private final boolean stableCrust;
	private final boolean subInterface;
	private final boolean subSlab;
	private final boolean volcanic;
	private final boolean grid;
	private final boolean faults;

	public NshmErf(Path path, Set<TectonicRegionType> trts,
			IncludeBackgroundOption gridOption) {
		this(HazardModel.load(path), trts, gridOption);
	}

	public NshmErf(HazardModel model, Set<TectonicRegionType> trts,
			IncludeBackgroundOption gridOption) {
		this.model = model;
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

	private void init() {

		// ERF initializers
		timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		timeSpan.setDuration(50.0);
		setTimeSpan(timeSpan);

		// nshmp-haz initializers
		Multimap<TectonicSetting, SourceTree> trees = model.trees();
		for (Entry<TectonicSetting, SourceTree> entry : trees.entries()) {

			TectonicSetting setting = entry.getKey();
			SourceTree tree = entry.getValue();
			SourceType type = tree.type();

			if (setting == TectonicSetting.SUBDUCTION) {
				if ((type == SourceType.INTERFACE || type == SourceType.INTERFACE_CLUSTER
						|| type == SourceType.INTERFACE_GRID || type == SourceType.INTERFACE_SYSTEM) && !subInterface) {
					continue;
				}
				if ((type == SourceType.SLAB || type == SourceType.INTRASLAB_GRID) && !subSlab) {
					continue;
				}
			}
			if (setting == TectonicSetting.STABLE_CRUST && !stableCrust) {
				continue;
			}
			if (setting == TectonicSetting.ACTIVE_CRUST && !activeCrust) {
				continue;
			}
			if (setting == TectonicSetting.VOLCANIC && !volcanic) {
				continue;
			}

			TectonicRegionType trt = NshmUtil.tectonicSettingToType(setting, type);
			List<NshmSource> sources = initTree(tree);
			sources.forEach(s -> s.setTectonicRegionType(trt));
			allSources.addAll(sources);
			sourceMap.putAll(trt, sources);
		}
	}

	public List<NshmSource> allSources() {
		return allSources;
	}

	public Multimap<TectonicRegionType, NshmSource> sourceMap() {
		return sourceMap;
	}

	private List<NshmSource> initTree(SourceTree tree) {
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

	private List<NshmSource> sourcesFromBranch(
			Branch<RuptureSet> branch,
			double duration) {

		RuptureSet ruptureSet = branch.value();
		double weight = branch.weight();

		switch (ruptureSet.type()) {

		case GRID:
			GriddedRuptureSet grs = (GriddedRuptureSet) ruptureSet;
			return (grid)
					? gridRuptureSetToSources(grs, weight, duration)
							: List.of();

		case ZONE:
			GriddedRuptureSet zrs = (GriddedRuptureSet) ruptureSet;
			return (faults)
					? gridRuptureSetToSources(zrs, weight, duration)
							: List.of();

		case INTERFACE_GRID:
			GriddedRuptureSet igrs = (GriddedRuptureSet) ruptureSet;
			return (subInterface && grid)
					? gridRuptureSetToSources(igrs, weight, duration)
							: List.of();

		case INTRASLAB_GRID:
			GriddedRuptureSet isgrs = (GriddedRuptureSet) ruptureSet;
			return (subSlab && grid)
					? gridRuptureSetToSources(isgrs, weight, duration)
							: List.of();

		case SLAB:
			GriddedRuptureSet slabRuptures = (GriddedRuptureSet) ruptureSet;
			return (subSlab && grid)
					? gridRuptureSetToSources(slabRuptures, weight, duration)
							: List.of();

		case FAULT_CLUSTER:
			ClusterRuptureSet crs = (ClusterRuptureSet) ruptureSet;
			return (faults)
					? clusterRuptureSetToSources(crs, weight, duration)
							: List.of();

		case FAULT_SYSTEM:
			SystemRuptureSet srs = (SystemRuptureSet) ruptureSet;
			return (faults)
					? systemRuptureSetToSources(srs, weight, duration)
							: List.of();

		case INTERFACE:
			return (subInterface && faults)
					? iterableRuptureSetToSources((IterableRuptureSet) ruptureSet, weight, duration)
							: List.of();

		case INTERFACE_CLUSTER:
			ClusterRuptureSet icrs = (ClusterRuptureSet) ruptureSet;
			return (subInterface && faults)
					? clusterRuptureSetToSources(icrs, weight, duration)
							: List.of();

		case INTERFACE_SYSTEM:
			SystemRuptureSet isrs = (SystemRuptureSet) ruptureSet;
			return (faults)
					? systemRuptureSetToSources(isrs, weight, duration)
							: List.of();

		default:
			return (faults)
					? iterableRuptureSetToSources((IterableRuptureSet) ruptureSet, weight, duration)
							: List.of();
		}
	}

	private static List<NshmSource> gridRuptureSetToSources(
			GriddedRuptureSet ruptureSet,
			double weight,
			double duration) {

		List<NshmSource> sources = new ArrayList<>();
		for (GridSource gridSource : ruptureSet) {
			sources.add(new NshmSource.Point(gridSource, weight, duration));
		}
		return sources;
	}

	private static List<NshmSource> iterableRuptureSetToSources(
			IterableRuptureSet ruptureSet,
			double weight,
			double duration) {

		return List.of(new NshmSource.Fault(ruptureSet, weight, duration));
	}

	private static List<NshmSource> systemRuptureSetToSources(
			SystemRuptureSet srs,
			double weight,
			double duration) {

		List<NshmSurface> surfaces = Arrays.stream(srs.sections)
				.map(section -> new NshmSurface(section))
				.collect(Collectors.toList());

		// SystemRuptureSet.stream() not supported but should be.
		// Iterator should work, even if it isn't used in nshm calc
		// pathways

		List<NshmSource> sources = new ArrayList<>(srs.size());
		for (int i = 0; i < srs.size(); i++) {
			SystemRupture source = (SystemRupture) srs.get(i);
			int[] sectionIndices = Indexing.bitsToIndices(source.bitset());
			List<NshmSurface> ruptureSurfaces = IntStream.of(sectionIndices)
					.mapToObj(surfaces::get)
					.collect(Collectors.toList());
			sources.add(new NshmSource.System(srs, source, weight, duration, ruptureSurfaces));
		}
		return sources;
	}

	private static List<NshmSource> clusterRuptureSetToSources(
			ClusterRuptureSet crs,
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
