package scratch.kevin.bbp;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;

import org.dom4j.Element;
import org.opensha.commons.metadata.XMLSaveable;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class BBP_Module implements XMLSaveable {
	
	public static final String VERSION = "17.3.0";
	
	public static enum VelocityModel {
		LA_BASIN_863("LA Basin 863 (m/s)", "LABasin863", "$BBP_INSTALL_GF/LABasin863/gp/genslip_nr_generic1d-gp01.vmod", 863),
		LA_BASIN_500("LA Basin 500 (m/s)", "LABasin500", "$BBP_INSTALL_GF/LABasin500/gp/nr02-vs500.fk1d", 500);
		
		private String name;
		private String xmlName;
		private String filePath;
		private double vs30;

		private VelocityModel(String name, String xmlName, String filePath, double vs30) {
			this.name = name;
			this.xmlName = xmlName;
			this.filePath = filePath;
			this.vs30 = vs30;
		}
		
		public String getDirName() {
			return xmlName;
		}
		
		public double getVs30() {
			return vs30;
		}
		
		public File getFilePath(String gfDir) {
			return new File(filePath.replaceAll(Matcher.quoteReplacement("$BBP_INSTALL_GF"), gfDir));
		}
		
		@Override
		public String toString() {
			return name;
		}
	}
	
	public static enum Method {
		GP("GP", 2d);
		
		private String xmlName;
		private double determLowpassFreq;

		private Method(String xmlName, double determLowpassFreq) {
			this.xmlName = xmlName;
			this.determLowpassFreq = determLowpassFreq;
		}
		
		public double getDetermLowpassFreq() {
			return determLowpassFreq;
		}
	}
	
	public static BBP_Module buildGenSlip(VelocityModel vm, File srcFile) {
		List<String> stagedFiles = Lists.newArrayList(vm.filePath, srcFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(fileName(vm.filePath)), stringArg(srcFile.getName()),
				stringArg(getGeneratedSRFName(srcFile)), stringArg(vm.xmlName));
		return new BBP_Module("Genslip", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildJBSim(VelocityModel vm, File srcFile, File srfFile, File siteFile) {
		String srfArg;
		List<String> stagedFiles;
		if (srfFile == null) {
			srfArg = getGeneratedSRFName(srcFile);
			stagedFiles = Lists.newArrayList(vm.filePath, srcFile.getAbsolutePath(), siteFile.getAbsolutePath());
		} else {
			srfArg = srfFile.getName();
			stagedFiles = Lists.newArrayList(vm.filePath, srcFile.getAbsolutePath(), srfFile.getAbsolutePath(),
					siteFile.getAbsolutePath());
		}
		List<Argument> arguments = Lists.newArrayList(stringArg(fileName(vm.filePath)), stringArg(srcFile.getName()),
				stringArg(srfArg), stringArg(siteFile.getName()), stringArg(vm.xmlName));
		return new BBP_Module("Jbsim", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildHFSims(VelocityModel vm, File srcFile, File srfFile, File siteFile) {
		String srfArg;
		List<String> stagedFiles;
		if (srfFile == null) {
			srfArg = getGeneratedSRFName(srcFile);
			stagedFiles = Lists.newArrayList(vm.filePath, srcFile.getAbsolutePath(), siteFile.getAbsolutePath());
		} else {
			srfArg = srfFile.getName();
			stagedFiles = Lists.newArrayList(vm.filePath, srcFile.getAbsolutePath(), srfFile.getAbsolutePath(),
					siteFile.getAbsolutePath());
		}
		List<Argument> arguments = Lists.newArrayList(stringArg(fileName(vm.filePath)), stringArg(srcFile.getName()),
				stringArg(srfArg), stringArg(siteFile.getName()), stringArg(vm.xmlName));
		return new BBP_Module("Hfsims", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildMatch(VelocityModel vm, File siteFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()), stringArg(vm.xmlName));
		List<KeywordArugment> kwArgs = Lists.newArrayList(new KeywordArugment("acc", "bool", "False"));
		return new BBP_Module("Match", stagedFiles, arguments, kwArgs);
	}
	
	public static BBP_Module buildPlotMap(File srcFile, File srfFile, File siteFile) {
		List<String> stagedFiles;
		List<Argument> arguments;
		if (srfFile == null) {
			stagedFiles = Lists.newArrayList(srcFile.getAbsolutePath(), siteFile.getAbsolutePath());
			arguments = Lists.newArrayList(stringArg(srcFile.getName()), stringArg(siteFile.getName()));
		} else {
			stagedFiles = Lists.newArrayList(srfFile.getAbsolutePath(), siteFile.getAbsolutePath());
			arguments = Lists.newArrayList(stringArg(srfFile.getName()), stringArg(siteFile.getName()));
		}
		return new BBP_Module("Plot_Map", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildPlotSeis(File siteFile, File srcFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath(), srcFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()), stringArg(srcFile.getName()),
				boolTrueArg, boolTrueArg);
		return new BBP_Module("PlotSeis", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildRotD50(File siteFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()));
		return new BBP_Module("RotD50", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildRotD100(File siteFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()));
		return new BBP_Module("RotD100", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildFAS(File siteFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()));
		return new BBP_Module("FAS", stagedFiles, arguments, null);
	}
	
	public static BBP_Module buildGenHTML(VelocityModel vm, Method method, File siteFile, File srcFile) {
		List<String> stagedFiles = Lists.newArrayList(siteFile.getAbsolutePath(), srcFile.getAbsolutePath());
		List<Argument> arguments = Lists.newArrayList(stringArg(siteFile.getName()), stringArg(srcFile.getName()),
				stringArg(vm.xmlName), new Argument("NoneType", "None"), stringArg(method.xmlName));
		return new BBP_Module("GenHTML", stagedFiles, arguments, null);
	}
	
	private static String fileName(String filePath) {
		if (filePath.contains("/"))
			return filePath.substring(filePath.lastIndexOf("/")+1);
		return filePath;
	}
	
	private static String getGeneratedSRFName(File srcFile) {
		String name = srcFile.getName();
		if (name.endsWith(".src"))
			name = name.substring(0, name.indexOf(".src"));
		return name+".srf";
	}
	
	private final String name;
	private final List<String> stagedFiles;
	private final List<Argument> arguments;
	private final List<KeywordArugment> kwArgs;
	
	private static Argument stringArg(String value) {
		return new Argument("str", value);
	}
	
	private static Argument boolTrueArg = new Argument("bool", "True");
	private static Argument boolFalseArg = new Argument("bool", "False");
	
	private static class Argument {
		final String type, value;
		
		public Argument(String type, String value) {
			this.type = type;
			this.value = value;
		}
	}
	
	private static class KeywordArugment {
		final String keyword, type, value;

		public KeywordArugment(String keyword, String type, String value) {
			this.keyword = keyword;
			this.type = type;
			this.value = value;
		}
	}
	
	private BBP_Module(String name, List<String> stagedFiles, List<Argument> arguments, List<KeywordArugment> kwArgs) {
		this.name = name;
		this.stagedFiles = Collections.unmodifiableList(stagedFiles);
		this.arguments = Collections.unmodifiableList(arguments);
		this.kwArgs = kwArgs;
	}

	@Override
	public Element toXMLMetadata(Element root) {
		Element modEl = root.addElement("BBP_Module");
		addTextEl(modEl, "name", name);
		Element stagedEl = modEl.addElement("staged_files");
		for (String file : stagedFiles)
			addTextEl(stagedEl, "file", file);
		Element argEl = modEl.addElement("arguments");
		for (Argument arg : arguments)
			addTextEl(argEl, "argument", arg.value, "type", arg.type);
		if (kwArgs != null) {
			Element kwEl = modEl.addElement("keyword_arguments");
			for (KeywordArugment arg : kwArgs)
				addTextEl(kwEl, "keyword_argument", arg.value, "keyword", arg.keyword, "type", arg.type);
		}
		return root;
	}
	
	private static Element addTextEl(Element root, String name, String text, String... attributePairs) {
		Element el = root.addElement(name);
		el.setText(text);
		if (attributePairs != null && attributePairs.length > 0) {
			Preconditions.checkState(attributePairs.length % 2 == 0);
			for (int i=0; i<attributePairs.length; i+=2)
				el.addAttribute(attributePairs[i], attributePairs[i+1]);
		}
		return el;
	}

}
