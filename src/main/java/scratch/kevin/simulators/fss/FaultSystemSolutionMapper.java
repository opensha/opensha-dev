package scratch.kevin.simulators.fss;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class FaultSystemSolutionMapper {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5892.instance();
		FaultSystemSolution compSol = catalog.getComparisonSolution();
		
		double overallMinMag = 6d;
		double sectFract = 0.5;
		int minSubSects = 2;
		
		double duration = catalog.getDurationYears();
		int skipYears;
		if (duration > 400000)
			skipYears = 20000;
		else if (duration > 100000)
			skipYears = 10000;
		else if (duration > 50000)
			skipYears = 5000;
		else
			skipYears = 2000;
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		mapper.setMinFractForInclusion(sectFract);
		List<RSQSimEvent> events = catalog.loader().skipYears(skipYears).minMag(overallMinMag)
				.minMappedSubSects(minSubSects, mapper).load();
		List<? extends FaultSection> subSects = catalog.getSubSects();
		
		System.out.println("Have "+events.size()+" events on "+subSects.size()+" sub sects");
		
		FaultSystemSolution sol = RSQSimUtils.buildFaultSystemSolution(subSects, catalog.getElements(), events, overallMinMag, mapper);
		
		if (compSol != null && compSol.hasModule(GridSourceProvider.class))
			sol.setGridSourceProvider(compSol.getGridSourceProvider());
		sol.write(new File(catalog.getCatalogDir(), "fss_m"+new DecimalFormat("0.#").format(overallMinMag)
				+"_skip"+skipYears+"_sectArea"+(float)sectFract+"_minSubSects"+minSubSects+".zip"));
	}

}
