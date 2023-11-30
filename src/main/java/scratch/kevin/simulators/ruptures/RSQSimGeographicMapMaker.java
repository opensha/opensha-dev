package scratch.kevin.simulators.ruptures;

import java.awt.BasicStroke;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

public class RSQSimGeographicMapMaker extends GeographicMapMaker {
	
	private RSQSimEvent event;
	private Color elemFillColor;
	private List<Double> elemScalars;
	private CPT elemCPT;
	private String elemScalarLabel;
	private boolean fillElemScalars;
	private Color elemOutlineColor;
	private float elemOutlineThickness;
	
	private Color eventHypoColor = null;
	private double eventHypoRadius = 0.2;
	

	public RSQSimGeographicMapMaker(Region region) {
		super(region);
		setWriteGeoJSON(false);
		setWritePDFs(false);
	}

	public RSQSimGeographicMapMaker(Region region, XY_DataSet[] politicalBoundaries) {
		super(region, politicalBoundaries);
		setWriteGeoJSON(false);
		setWritePDFs(false);
	}
	
	public void plotEvent(RSQSimEvent event, Color fillColor, Color outlineColor, float outlineThickness) {
		clearEvent();
		Preconditions.checkState(fillColor != null || outlineColor != null);
		this.event = event;
		this.elemFillColor = fillColor;
		this.elemOutlineColor = outlineColor;
		this.elemOutlineThickness = outlineThickness;
	}
	
	public void plotEventFillScalars(RSQSimEvent event, List<Double> scalars, CPT cpt, String label) {
		plotEventFillScalars(event, scalars, cpt, null, Float.NaN, label);
	}
	
	public void plotEventFillScalars(RSQSimEvent event, List<Double> scalars, CPT cpt, Color outlineColor, float outlineThickness, String label) {
		clearEvent();
		Preconditions.checkState(scalars.size() == event.getNumElements());
		Preconditions.checkNotNull(cpt);
		this.event = event;
		this.elemScalars = scalars;
		this.elemCPT = cpt;
		this.elemScalarLabel = label;
		this.elemOutlineColor = outlineColor;
		this.elemOutlineThickness = outlineThickness;
		this.fillElemScalars = true;
	}
	
	public void plotEventOutlineScalars(RSQSimEvent event, List<Double> scalars, CPT cpt, float outlineThickness, String label) {
		clearEvent();
		Preconditions.checkState(scalars.size() == event.getNumElements());
		Preconditions.checkNotNull(cpt);
		this.event = event;
		this.elemScalars = scalars;
		this.elemScalarLabel = label;
		this.elemCPT = cpt;
		this.elemOutlineThickness = outlineThickness;
		this.fillElemScalars = false;
	}
	
	private void clearEvent() {
		this.event = null;
		this.elemFillColor = null;
		this.elemScalars = null;
		this.elemCPT = null;
		this.elemScalarLabel = null;
		this.elemOutlineColor = null;
		this.elemOutlineThickness = Float.NaN;
	}
	
	public void plotEventHypocenter(Color color) {
		plotEventHypocenter(color, eventHypoRadius);
	}
	
	public void plotEventHypocenter(Color color, double radius) {
		eventHypoColor = color;
		eventHypoRadius = radius;
	}

	@Override
	protected PlotBuilder initPlotBuilder() {
		return new RSQSimPlotBuilder();
	}
	
	protected class RSQSimPlotBuilder extends PlotBuilder {

		@Override
		protected void plotFirst() {
			super.plotFirst();
			
			if (event != null) {
				// plot event
				
				ArrayList<SimulatorElement> elems = event.getAllElements();
				
				XY_DataSet[] elemXYs = buildElemXYs(elems);
				
				if (elemCPT != null) {
					Preconditions.checkState(elems.size() == elemScalars.size());
					
					if (fillElemScalars) {
						// fill them
						for (int i=0; i<elemXYs.length; i++) {
							double val = elemScalars.get(i);
							if (skipNaNs && Double.isNaN(val))
								continue;
							funcs.add(elemXYs[i]);
							chars.add(new PlotCurveCharacterstics(
									PlotLineType.POLYGON_SOLID, 1f, elemCPT.getColor((float)val)));
						}
						if (elemOutlineColor != null) {
							// add outlines on top
							PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(
									PlotLineType.SOLID, elemOutlineThickness, elemOutlineColor);
							for (int i=0; i<elemXYs.length; i++) {
								funcs.add(elemXYs[i]);
								chars.add(outlineChar);
							}
						}
					} else {
						// outlines
						for (int i=0; i<elemXYs.length; i++) {
							double val = elemScalars.get(i);
							if (skipNaNs && Double.isNaN(val))
								continue;
							funcs.add(elemXYs[i]);
							chars.add(new PlotCurveCharacterstics(
									PlotLineType.SOLID, elemOutlineThickness, elemCPT.getColor((float)val)));
						}
					}
					if (elemScalarLabel != null)
						cptLegend.add(buildCPTLegend(elemCPT, elemScalarLabel));
				} else {
					// constant color
					if (elemFillColor != null) {
						PlotCurveCharacterstics fillChar = new PlotCurveCharacterstics(
								PlotLineType.POLYGON_SOLID, 1f, elemFillColor);
						for (int i=0; i<elemXYs.length; i++) {
							funcs.add(elemXYs[i]);
							chars.add(fillChar);
						}
					}
					if (elemOutlineColor != null) {
						PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(
								PlotLineType.SOLID, elemOutlineThickness, elemOutlineColor);
						for (int i=0; i<elemXYs.length; i++) {
							funcs.add(elemXYs[i]);
							chars.add(outlineChar);
						}
					}
				}
			}
			
			if (eventHypoColor != null) {
				double outerRatio = 2.618;
				int num = 10;
				double radsEach = Math.PI/5d;
				
				Location loc = RSQSimUtils.getHypocenter(event);
				
				DefaultXY_DataSet poly = new DefaultXY_DataSet();
				
				Range xRange = getXRange();
				Range yRange = getYRange();
				double plotAspectRatio = PlotUtils.calcAspectRatio(xRange, yRange, true);
				double aspect = plotAspectRatio*xRange.getLength()/yRange.getLength();
//				System.out.println("plotAspectRatio="+plotAspectRatio);
//				System.out.println("aspect="+aspect);
				
				for (int i=0; i<10; i++) {
					double dist;
					if (i % 2 == 1)
						dist = eventHypoRadius;
					else
						dist = eventHypoRadius/outerRatio;
					double angle = radsEach*i;
					
					double dx = Math.sin(angle)*dist/aspect;
					double dy = Math.cos(angle)*dist;
					
					double x = loc.lon + dx;
					double y = loc.lat + dy;
					poly.set(x, y);
				}
				poly.set(poly.get(0));
				
				funcs.add(poly);
				chars.add(new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, eventHypoColor));
				funcs.add(poly);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
			}
		}
		
	}
	
	private static XY_DataSet[] buildElemXYs(List<? extends SimulatorElement> elems) {
		XY_DataSet[] ret = new XY_DataSet[elems.size()];
		
		for (int i=0; i<ret.length; i++) {
			SimulatorElement elem = elems.get(i);
			Vertex[] vertexes = elem.getVertices();
			double[] xVals = new double[vertexes.length+1];
			double[] yVals = new double[vertexes.length+1];
			for (int n=0; n<xVals.length; n++) {
				xVals[n] = vertexes[n%vertexes.length].getLongitude();
				yVals[n] = vertexes[n%vertexes.length].getLatitude();
			}
			ret[i] = new DefaultXY_DataSet(xVals, yVals);
		}
		
		return ret;
	}

}
