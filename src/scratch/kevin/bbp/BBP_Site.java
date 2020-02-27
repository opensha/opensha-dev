package scratch.kevin.bbp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class BBP_Site {
	
	private String name;
	private Location loc;
	private double vs30;
	private double loPassFreq;
	private double hiPassFreq;

	public BBP_Site(String name, Location loc, double vs30, double loPassFreq, double hiPassFreq) {
		this.name = name;
		this.loc = loc;
		this.vs30 = vs30;
		this.loPassFreq = loPassFreq;
		this.hiPassFreq = hiPassFreq;
	}

	public static void writeToFile(File file, Collection<BBP_Site> sites) throws IOException {
		FileWriter fw = new FileWriter(file);
		
		fw.write("#SLong    SLat     RSN   Vs30(m/s) LoPass_Freq(Hz) HiPass_Freq(Hz)\n");
		for (BBP_Site site : sites)
			fw.write((float)site.loc.getLongitude()+" "+(float)site.loc.getLatitude()+" "+site.name+" "
					+(float)site.vs30+" "+(float)site.loPassFreq+" "+(float)site.hiPassFreq+"\n");
		
		fw.close();
	}
	
	public static List<BBP_Site> readFile(File file) throws IOException {
		if (file.isDirectory()) {
			for (File f : file.listFiles())
				if (f.getName().toLowerCase().endsWith(".stl"))
					return BBP_Site.readFile(f);
			throw new FileNotFoundException("No .stl files found in "+file.getAbsolutePath());
		}
		List<BBP_Site> sites = new ArrayList<>();
		
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			String[] split = line.split("\\s+");
			Preconditions.checkState(split.length == 6);
			Location loc = new Location(Double.parseDouble(split[1]), Double.parseDouble(split[0]));
			String name = split[2];
			double vs30 = Double.parseDouble(split[3]);
			double loPassFreq = Double.parseDouble(split[4]);
			double hiPassFreq = Double.parseDouble(split[5]);
			
			sites.add(new BBP_Site(name, loc, vs30, loPassFreq, hiPassFreq));
		}
		
		return sites;
	}

	public String getName() {
		return name;
	}

	public Location getLoc() {
		return loc;
	}

	public double getVs30() {
		return vs30;
	}

	public double getLoPassFreq() {
		return loPassFreq;
	}

	public double getHiPassFreq() {
		return hiPassFreq;
	}
	
	public Site buildGMPE_Site(VelocityModel vm) {
		Site gmpeSite = new Site(getLoc(), getName());
		
		Vs30_Param vs30Param = new Vs30_Param(vs30);
		vs30Param.setValueAsDefault();
		gmpeSite.addParameter(vs30Param);
		
		Vs30_TypeParam vs30TypeParam = new Vs30_TypeParam();
		vs30TypeParam.setValue(Vs30_TypeParam.VS30_TYPE_INFERRED); // TODO
		gmpeSite.addParameter(vs30TypeParam);
		
		DepthTo1pt0kmPerSecParam z10Param = new DepthTo1pt0kmPerSecParam();
		z10Param.setValue(vm == null ? null : vm.getZ10()*1000d);
		gmpeSite.addParameter(z10Param);
		
		DepthTo2pt5kmPerSecParam z25Param = new DepthTo2pt5kmPerSecParam();
		z25Param.setValue(vm == null ? null : vm.getZ25());
		gmpeSite.addParameter(z25Param);
		
		return gmpeSite;
	}

}
