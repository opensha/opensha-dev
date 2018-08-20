package scratch.aftershockStatistics.aafs;

import java.util.List;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;

import scratch.aftershockStatistics.ComcatException;
import scratch.aftershockStatistics.CompactEqkRupList;

import scratch.aftershockStatistics.pdl.PDLProductBuilderOaf;
import scratch.aftershockStatistics.pdl.PDLSender;
import gov.usgs.earthquake.product.Product;


/**
 * Support functions for PDL reporting.
 * Author: Michael Barall 06/24/2018.
 */
public class PDLSupport extends ServerComponent {




	//----- Task execution subroutines : PDL operations -----




	// Send a report to PDL.
	// Throw an exception if the report failed.

	public void send_pdl_report (TimelineStatus tstatus) throws Exception {

		// The JSON file to send

		String jsonText = tstatus.forecast_results.get_pdl_model();

		// The event network and code

		String eventNetwork = tstatus.forecast_mainshock.mainshock_network;
		String eventCode = tstatus.forecast_mainshock.mainshock_code;

		// The event ID, which for us identifies the timeline

		String eventID = sg.alias_sup.timeline_id_to_pdl_code (tstatus.event_id);

		// Modification time, 0 means now

		long modifiedTime = 0L;

		// Review status, false means automatically generated

		boolean isReviewed = false;

		// Build the product

		Product product = PDLProductBuilderOaf.createProduct (eventID, eventNetwork, eventCode, isReviewed, jsonText, modifiedTime);

		// Sign the product

		PDLSender.signProduct(product);

		// Send the product, true means it is text

		PDLSender.sendProduct(product, true);

		return;
	}




	// Return true if this machine is primary for sending reports to PDL, false if secondary

	public boolean is_pdl_primary () {

		// For now, just assume primary

		return true;
	}



	//----- Construction -----


	// Default constructor.

	public PDLSupport () {}

}
