package scratch.kevin.ucerf3.downDipSubSectTest;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.dom4j.Element;
import org.opensha.commons.util.IDPairing;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.inversion.laughTest.AbstractPlausibilityFilter;
import scratch.UCERF3.inversion.laughTest.PlausibilityConfiguration;
import scratch.UCERF3.inversion.laughTest.UCERF3PlausibilityConfig;

public class DownDipTestPlausibilityConfig implements PlausibilityConfiguration {
	
	private DownDipSubSectBuilder downDipBuilder;
	private UCERF3PlausibilityConfig u3Config;
	
	private List<AbstractPlausibilityFilter> filters;

	public DownDipTestPlausibilityConfig(DownDipSubSectBuilder downDipBuilder) {
		this.downDipBuilder = downDipBuilder;
		u3Config = UCERF3PlausibilityConfig.getDefault();
		u3Config.setCoulombFilter(null); // no coulomb
		// disable all of the azimuthal filters as they don't work with down-dip subsections
		u3Config.setMaxCmlAzimuthChange(Double.NaN);
		u3Config.setMaxTotAzimuthChange(Double.NaN);
		u3Config.setMaxAzimuthChange(Double.NaN);
	}

	@Override
	public Element toXMLMetadata(Element root) {
		// TODO: this currently doesn't store anything to recreate this test configuration,
		// but that shouldn't matter (since it's just a test)
		return u3Config.toXMLMetadata(root);
	}

	@Override
	public synchronized List<AbstractPlausibilityFilter> buildPlausibilityFilters(Map<IDPairing, Double> azimuths,
			Map<IDPairing, Double> distances, List<List<Integer>> sectionConnectionsListList,
			List<? extends FaultSection> subSectData) {
		filters = new ArrayList<>();
		filters.addAll(u3Config.buildPlausibilityFilters(azimuths, distances,
				sectionConnectionsListList, subSectData));
		filters.add(new RectangularityFilter(downDipBuilder, 2));
		return filters;
	}

	@Override
	public List<AbstractPlausibilityFilter> getPlausibilityFilters() {
		return filters;
	}

	@Override
	public double getMaxJumpDist() {
		return u3Config.getMaxJumpDist();
	}

}
