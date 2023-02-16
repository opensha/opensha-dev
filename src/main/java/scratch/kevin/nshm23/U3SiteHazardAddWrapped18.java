package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SiteLogicTreeHazardCurveCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.NEHRP_TestCity;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.NshmErf;

public class U3SiteHazardAddWrapped18 {

	public static void main(String[] args) throws IOException {
		if (args.length != 4) {
			System.err.println("USAGE: <source-dir> <output-zip> <model-path> <INCLUDE/EXCLUDE/ONLY>");
			System.exit(1);
		}
		
		File sourceDir = new File(args[0]);
		File outputZip = new File(args[1]);
		Path modelPath = Path.of(args[2]);
		IncludeBackgroundOption griddedOp = IncludeBackgroundOption.valueOf(args[3]);
		
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW, TectonicRegionType.STABLE_SHALLOW);
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}

		NshmErf erf = new NshmErf(modelPath, trts, griddedOp);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.get();
		gmpe.setParamDefaults();
		
		double[] periods = { 0d, 1d };
		
		DiscretizedFunc[] perXVals = new DiscretizedFunc[periods.length];
		DiscretizedFunc[] perLogXVals = new DiscretizedFunc[periods.length];
		
		for (int p=0; p<periods.length; p++) {
			DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(periods[p] > 0 ? SA_Param.NAME : PGA_Param.NAME);
			DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : xVals)
				logXVals.set(Math.log(pt.getX()), 0d);
			perXVals[p] = xVals;
			perLogXVals[p] = logXVals;
		}
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		calc.setMaxSourceDistance(500d);
		
		ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(outputZip));
		
		CSVFile<String> sitesCSV = new CSVFile<>(true);
		sitesCSV.addLine("Name", "Latitude", "Longitude");
		
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		for (NEHRP_TestCity city : NEHRP_TestCity.values()) {
			Location loc = city.location();
			if (!reg.contains(loc))
				continue;
			Site site = new Site(loc);
			String siteName = city.toString();
			site.addParameterList(gmpe.getSiteParams());

			System.out.println("Calculating for "+city);
			
			sitesCSV.addLine(siteName, (float)loc.getLatitude()+"", (float)loc.getLongitude()+"");
			
			for (int p=0; p<periods.length; p++) {
				String csvName = siteName.replaceAll("\\W+", "_");
				if (periods[p] == 0d) {
					gmpe.setIntensityMeasure(PGA_Param.NAME);
					csvName += "_pga.csv";
				} else {
					gmpe.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), periods[p]);
					csvName += "_sa_"+(float)periods[p]+".csv";
				}
				
				File srcFile = new File(sourceDir, csvName);
				CSVFile<String> csv;
				if (srcFile.exists()) {
					System.out.println("Using UCERF3 file: "+csvName);
					csv = CSVFile.readFile(srcFile, true);
				} else {
					DiscretizedFunc logXVals = perLogXVals[p];
					DiscretizedFunc xVals = perXVals[p];
					calc.getHazardCurve(logXVals, site, gmpe, erf);

					DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
					for (int i=0; i<xVals.size(); i++)
						curve.set(xVals.getX(i), logXVals.getY(i));

					System.out.println("Hazard curve:\n"+curve);
					System.out.println("2 in 50: "+curve.getFirstInterpolatedX_inLogXLogYDomain(ReturnPeriods.TWO_IN_50.oneYearProb));
					
					List<String> header = new ArrayList<>();
					header.add("Site Name");
					header.add("Branch Index");
					header.add("Branch Weight");
					for (Point2D pt : xVals)
						header.add((float)pt.getX()+"");
					
					csv = new CSVFile<>(true);
					csv.addLine(header);
					List<String> line = new ArrayList<>();
					line.add(siteName);
					line.add("0");
					line.add("1.0");
					for (Point2D pt : curve)
						line.add(pt.getY()+"");
					csv.addLine(line);
				}
				
				zout.putNextEntry(new ZipEntry(csvName));
				csv.writeToStream(zout);
				zout.flush();
				zout.closeEntry();
			}
		}
		// write sites
		zout.putNextEntry(new ZipEntry(MPJ_SiteLogicTreeHazardCurveCalc.SITES_CSV_FILE_NAME));
		sitesCSV.writeToStream(zout);
		zout.flush();
		zout.closeEntry();
		
		zout.close();
		
	}

}
