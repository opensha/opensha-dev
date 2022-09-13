package gov.usgs.earthquake.nshmp.model;

import static java.util.stream.Collectors.toList;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.opensha.commons.data.TimeSpan;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;

import gov.usgs.earthquake.nshmp.Maths;
import gov.usgs.earthquake.nshmp.data.Indexing;
import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemSource;
import gov.usgs.earthquake.nshmp.tree.Branch;

public class NshmErf extends AbstractERF {

  private final HazardModel model;
  private final List<NshmSource> allSources;
  private final Multimap<TectonicRegionType, NshmSource> sourceMap;

  public NshmErf(Path path) {
    model = HazardModel.load(path);
    allSources = new ArrayList<>();
    sourceMap = MultimapBuilder
        .enumKeys(TectonicRegionType.class)
        .arrayListValues()
        .build();
    init();
    System.out.println(allSources.size());
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
    return sources;
  }

  private static List<NshmSource> sourcesFromBranch(
      Branch<RuptureSet<? extends Source>> branch,
      double duration) {

    RuptureSet<? extends Source> ruptureSet = branch.value();
    double weight = branch.weight();

    switch (ruptureSet.type()) {

      // custom cluster and system handlers
      case FAULT_CLUSTER:
        ClusterRuptureSet crs = (ClusterRuptureSet) ruptureSet;
        return clusterRuptureSetToSources(crs, weight, duration);
      case FAULT_SYSTEM:
        SystemRuptureSet srs = (SystemRuptureSet) ruptureSet;
        return systemRuptureSetToSources(srs, weight, duration);

      default:
        return ruptureSetToSources(ruptureSet, weight, duration);
    }
  }

  private static List<NshmSource> ruptureSetToSources(
      RuptureSet<? extends Source> ruptureSet,
      double weight,
      double duration) {

    return ruptureSet.stream()
        .map(src -> new NshmSource(src, weight, duration))
        .collect(toList());
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
      SystemSource source = srs.get(i);
      int[] sectionIndices = Indexing.bitsToIndices(source.bitset());
      List<NshmSurface> ruptureSurfaces = IntStream.of(sectionIndices)
          .mapToObj(surfaces::get)
          .collect(Collectors.toList());
      sources.add(new NshmSource(source, ruptureSurfaces, weight, duration));
    }
    return sources;
  }

  private static List<NshmSource> clusterRuptureSetToSources(
      ClusterRuptureSet crs,
      double weight,
      double duration) {

    ClusterSource cs = crs.get(0); // known to be only one
    double rate = cs.rate();
    return cs.ruptureSets().stream()
        .map(rs -> (FaultRuptureSet) rs)
        .map(frs -> new NshmSource(frs.get(0), weight * rate, duration))
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
    allSources.stream()
        .map(NshmSource::ruptures)
        .flatMap(List::stream)
        .forEach(rup -> rup.setProbability(
            Maths.rateToProbability(rup.rate * rup.weight, duration)));
  }

  @Override
  public String getName() {
    return model.name();
  }

}
