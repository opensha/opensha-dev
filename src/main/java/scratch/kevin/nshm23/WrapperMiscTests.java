package scratch.kevin.nshm23;

import java.io.File;
import java.util.Collection;
import java.util.Set;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.util.TectonicRegionType;

import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.nshmp.shaded.model.NshmErf;
import org.opensha.nshmp.shaded.model.NshmpSourceTree;
import org.opensha.nshmp.shaded.model.NshmpTectonicSetting;

public class WrapperMiscTests {

	public static void main(String[] args) {
		File modelsDir = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models");
		
		File modelDir = new File(modelsDir, "nshm-conus-6.0.0");
		
		NshmpHazardModel model = NshmpHazardModel.load(modelDir.toPath());
//		NshmpTectonicSetting setting = NshmpTectonicSetting.ACTIVE_CRUST;
//		Collection<NshmpSourceTree> trees = model.trees().get(setting);
//		System.out.println("Iterating over "+trees.size()+" trees");
//		for (NshmpSourceTree tree : trees) {
//			System.out.println(tree);
//			tree.
//		}
		
		NshmErf erf = new NshmErf(model, Set.of(TectonicRegionType.ACTIVE_SHALLOW), IncludeBackgroundOption.ONLY);
		erf.updateForecast();
		
		System.out.println("have "+erf.getNumSources()+" sources");
		for (ProbEqkSource source : erf) {
			System.out.println("Source "+source.getName()+" of type "+source.getClass().getName()+" with "+source.getNumRuptures()+" rups");
			for (ProbEqkRupture rup : source) {
				RuptureSurface surf = rup.getRuptureSurface();
				System.out.println("M"+rup.getMag()+", rake="+rup.getAveRake()+", dip="+surf.getAveDip()
					+", surfType="+surf.getClass().getName());
				Location loc0 = surf.getEvenlyDiscretizedLocation(0);
				Location offset = LocationUtils.location(loc0, 0d, 10d);
				System.out.println("\trEpi=10: rJB="+(float)surf.getDistanceJB(offset)+", rRup="+(float)surf.getDistanceRup(offset)+", rX="+(float)surf.getDistanceX(offset));
			}
			break;
		}
	}

}
