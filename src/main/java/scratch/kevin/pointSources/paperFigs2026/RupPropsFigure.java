package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.PointSource.FocalMechRuptureSurfaceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.util.FocalMech;

import net.mahdilamb.colormap.Colors;

public class RupPropsFigure {

	public static void main(String[] args) throws IOException {
		Range magRange = new Range(3.5, 7.9);
		
		EvenlyDiscretizedFunc magBins = new EvenlyDiscretizedFunc(magRange.getLowerBound(), magRange.getUpperBound(), 1000);
		
		FocalMechRuptureSurfaceBuilder origBuilder = PointSourceNshm.SURF_BUILDER_DEFAULT;
		FocalMechRuptureSurfaceBuilder modBuilder = UPDATED_SURF_BUILDER;
		
		EvenlyDiscretizedFunc origZTorFunc = magBins.deepClone();
		EvenlyDiscretizedFunc origZBotDipFunc = magBins.deepClone();
		EvenlyDiscretizedFunc origZBotSSFunc = magBins.deepClone();
		EvenlyDiscretizedFunc origLenFunc = magBins.deepClone();
		
		EvenlyDiscretizedFunc modZTorDipFunc = magBins.deepClone();
		EvenlyDiscretizedFunc modZTorSSFunc = magBins.deepClone();
		EvenlyDiscretizedFunc modZBotDipFunc = magBins.deepClone();
		EvenlyDiscretizedFunc modZBotSSFunc = magBins.deepClone();
		EvenlyDiscretizedFunc modLenDipFunc = magBins.deepClone();
		EvenlyDiscretizedFunc modLenSSFunc = magBins.deepClone();
		Location loc = new Location(0, 0);
		for (int i=0; i<magBins.size(); i++) {
			double mag = magBins.getX(i);
			RuptureSurface origSurfSS = origBuilder.getSurface(loc, mag, FocalMech.STRIKE_SLIP, 0);
			RuptureSurface origSurfDip = origBuilder.getSurface(loc, mag, FocalMech.REVERSE, 0);
			RuptureSurface modSurfSS = modBuilder.getSurface(loc, mag, FocalMech.STRIKE_SLIP, 0);
			RuptureSurface modSurfDip = modBuilder.getSurface(loc, mag, FocalMech.REVERSE, 0);
			
			origZTorFunc.set(i, origSurfSS.getAveRupTopDepth());
			modZTorDipFunc.set(i, modSurfDip.getAveRupTopDepth());
			modZTorSSFunc.set(i, modSurfSS.getAveRupTopDepth());
			
			origZBotDipFunc.set(i, origSurfDip.getAveRupBottomDepth());
			origZBotSSFunc.set(i, origSurfSS.getAveRupBottomDepth());
			modZBotDipFunc.set(i, modSurfDip.getAveRupBottomDepth());
			modZBotSSFunc.set(i, modSurfSS.getAveRupBottomDepth());
			
			origLenFunc.set(i, origSurfSS.getAveLength());
			modLenDipFunc.set(i, modSurfDip.getAveLength());
			modLenSSFunc.set(i, modSurfSS.getAveLength());
		}
		
		Color origColor = Colors.tab_orange;
		Color modColor = Colors.tab_blue;
		
		List<DiscretizedFunc> depthFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> depthChars = new ArrayList<>();
		
		origZTorFunc.setName("NSHM23 as-published (strike-slip)");
		depthFuncs.add(origZTorFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, origColor));
		origZBotSSFunc.setName("NSHM23 as-published (dipping)");
		depthFuncs.add(origZBotSSFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, origColor));
		depthFuncs.add(origZBotDipFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, origColor));
		
		modZTorSSFunc.setName("Smooth Ztor transition (strike-slip)");
		depthFuncs.add(modZTorSSFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, modColor));
		modZTorDipFunc.setName("Smooth Ztor transition (dipping)");
		depthFuncs.add(modZTorDipFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, modColor));
		depthFuncs.add(modZBotSSFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, modColor));
		depthFuncs.add(modZBotDipFunc);
		depthChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, modColor));
		
		PlotSpec depthPlot = new PlotSpec(depthFuncs, depthChars, " ", "Magnitude", "Depth (km)");
		depthPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		depthPlot.setYAxisInverted(true);
		
		List<DiscretizedFunc> lengthFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> lengthChars = new ArrayList<>();
		
		origLenFunc.setName("NSHM23 as-published (WC94 SRL)");
		lengthFuncs.add(origLenFunc);
		lengthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, origColor));
		
		modLenSSFunc.setName("Leonard (2010) strike-slip");
		lengthFuncs.add(modLenSSFunc);
		lengthChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, modColor));
		
		modLenDipFunc.setName("Leonard (2010) dipping");
		lengthFuncs.add(modLenDipFunc);
		lengthChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, modColor));
		
		PlotSpec lengthPlot = new PlotSpec(lengthFuncs, lengthChars, " ", "Magnitude", "Length (km)");
		lengthPlot.setLegendInset(RectangleAnchor.TOP_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(List.of(depthPlot, lengthPlot), List.of(false), List.of(false, false), List.of(magRange), null);
		
		PlotUtils.writePlots(FIGURES_DIR, "rupture_properties", gp, 1000, 700, true, true, false);
	}

}
