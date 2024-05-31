package scratch.kevin.nshm23.prvi.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.faultSurface.FaultSection;

public class SourceZoneMap {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		
		PRVI25_CrustalFaultModels fm = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1;
		PRVI25_CrustalDeformationModels dm = PRVI25_CrustalDeformationModels.GEOLOGIC;
		List<? extends FaultSection> subSects = dm.build(fm);
		subSects = new ArrayList<>(subSects);
		for (int i=subSects.size(); --i>=0;)
			if (!subSects.get(i).isProxyFault())
				subSects.remove(i);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(subSects);
		mapMaker.setSectPolygonChar(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, new Color(0, 0, 0, 60)));
		mapMaker.plot(outputDir, "proxy_zones", "PRVI Source Zones and Proxy Faults");
	}

}
