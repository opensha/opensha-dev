package scratch.aftershockStatistics.pdl;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Map;

import org.json.simple.JSONObject;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.distribution.SocketProductSender;
import gov.usgs.earthquake.product.ByteContent;
import gov.usgs.earthquake.product.Content;
import gov.usgs.earthquake.product.Product;
import gov.usgs.earthquake.product.ProductId;
import scratch.aftershockStatistics.USGS_AftershockForecast;

public class OAF_Publisher {
	
	public static Product createProduct(String eventID, USGS_AftershockForecast model) {
		String source = "opensha-oaf-"+eventID;
		String type = "oaf";
		String code = source+":"+type;
		
		Product product;
		ProductId productId;

		// Need an ID for the product. The ID is for the product only and does
		// not necessarily reflect "event" properties, however, in practice
		// they are typically quite similar

		// A product ID uniquely identifies a product based on "source", "type",
		// and "code" properties. An "updateTime" property indicates a specific
		// version of a product based on when the product was created. By default,
		// newer versions of the same product (source/type/code) will replace
		// older versions.

		// Read more about source, type, code here:
		// http://ehppdl1.cr.usgs.gov/userguide/products/index.html
		productId = new ProductId(
				source,   // source
				type,     // type
				code      // code
				// updateTime - Defaults to NOW if not specificied, typically, do not
				//              specify this as a best practice
				);

		// Use Product.STATUS_UPDATE to create or update a product.
		// Use Product.STATUS_DELETE to remove this product.
		product = new Product(
				productId,
				Product.STATUS_UPDATE
				);

		// This is a legacy vestige that does not get used anymore, but we need
		// to add it otherwise sending will fail. Anything will do ...
		try {
			product.setTrackerURL(new URL("http://www.google.com/"));
		} catch (MalformedURLException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}


		attachPropertiesToProduct(product, eventID); // simply key-value properties
		attachContentToProduct(product, model);    // file- or byte-based content

		return product;
	}

	/**
	 * Attaches properties (key-value pairs) to a product.
	 *
	 * @param product
	 *     The product on to which the properties should be attached.
	 */
	private static void attachPropertiesToProduct (Product product, String eventID) {
		Map<String, String> properties;

		// Properties is a reference to the product's internal properties hash,
		// so updates made here are immediately reflected on the product properties
		properties = product.getProperties();


		// Typically, you always want to set the eventsource and eventsourcecode
		// properties. These are the strongest way to ensure the product being
		// sent will associate to the correct event once it reaches ComCat
		
		// detect event source
		Preconditions.checkState(!eventID.isEmpty());
		int numCharsBeforeDigit = 0;
		boolean hasDigit = false;
		for (int i=0; i<eventID.length(); i++) {
			if (Character.isDigit(eventID.charAt(i))) {
				hasDigit = true;
				break;
			}
			numCharsBeforeDigit++;
		}
		Preconditions.checkState(hasDigit && numCharsBeforeDigit>0, "couldn't detect event source from: %s", eventID);
		
		String eventSource = eventID.substring(0, numCharsBeforeDigit);
		String eventCode = eventID.substring(numCharsBeforeDigit);
		
		System.out.println("eventSource='"+eventSource+"', eventCode='"+eventCode+"'");

		// Note: We could also use the NC information "nc" and "72689331"
		// or any other eventid that references the same event that actually
		// occurred...
		properties.put(Product.EVENTSOURCE_PROPERTY, eventSource);
		properties.put(Product.EVENTSOURCECODE_PROPERTY, eventCode);


		// You may also want to set these other event-related properties
		// which are used as a fall-back for association purposes. There are pros
		// and cons to doing so that I won't discuss here. Dealer's choice!

		// properties.put(Product.EVENTTIME_PROPERTY, "2016-09-03T03:27:57.170UTC");
//		 properties.put(Product.MAGNITUDE_PROPERTY, "5.7");
		// properties.put(Product.LATITUDE_PROPERTY, "40.406");
		// properties.put(Product.LONGITUDE_PROPERTY, "-125.469");

		// properties.put(Product.DEPTH_PROPERTY, "5.6");

		// Note: Product.VERSION_PROPERTY does not indicate a "version" within
		// ComCat, but is really just an identifier for senders/network operators.
		// Typically it is not recommended to send this property.

		// properties.put(Proudct.VERSION_PROPERTY, "01");


		// This property is interpreted on the event page and dicates whether
		// an orange "X" or green checkmark is displayed at the top of the
		// product page. This status icon informs the user if this is an
		// automatically-generated forecast, or if a scientist manually generated
		// it. If this property is omitted, it is assumed the product was
		// automatically generated.
		properties.put("review-status", "reviewed"); // "reviewed" or "automatic"


		// Also add properties specific to your product-type, in this case,
		// things that are specific to the aftershock forecast.
		// We may want to discuss what specifically you add here...

//		properties.put("wcradius", "167.55"); // TODO
		// properties.put("key", "value"); // ... as needed/appropriate ...
	}
	
	/**
	 * If you have data structures/objects in runtime memory and want to attach
	 * them to the product, this is a good way to do that.
	 *
	 * @param product
	 *     The product on to which the content should be attached.
	 */
	private static void attachContentToProduct (Product product, USGS_AftershockForecast model) {
		Map<String, Content> contents = product.getContents();

//		// This content will be stored on disk in ComCat. For the event pages
//		// to access this content, they must make a separate HTTP request (which
//		// can slow page rendering). The only thing that makes this "disk"
//		// content is that it receives a non-empty content-path when
//		// it is added to the contents map. It has nothing to do with the
//		// origination of the content itself... It is perfectly fine for this
//		// content to be binary blobs if necessary.
//		diskContent = new ByteContent("Here is some disk content".getBytes());
//		diskContent.setContentType(/*String*/ null); // `null` implies "text/plain"
//		diskContent.setLastModified(/*Date*/ null); // `null` implies NOW
//		contents.put("disk-content.txt", diskContent);

		// This content will be embedded into the event details feed from ComCat
		// and is immediately available to the event page. The only thing that
		// makes this "inline" content is that the content-path is an empty string.
		// Content added as inline should be text-based and _not_ binary blobs.
		JSONObject json = model.buildJSON();
		ByteContent inlineContent = new ByteContent(
//				"{\"key\": \"Here is some inline content\"}".getBytes());
				json.toJSONString().getBytes());
		inlineContent.setContentType("application/json");
		inlineContent.setLastModified(/*Date*/ null); // `null` implies now
		contents.put("", inlineContent); // <-- Inline content because empty string

		// In general, you can serialize any Java object to a String and then
		// use String.getBytes() to generate ByteContent for a product.
	}
	
	public static void sendProduct(Product product) throws Exception {
		SocketProductSender sender;

		// SocketProductSenders send directly to a PDL HUB and do not introduce
		// any polling latency.
		sender = new SocketProductSender("pdldevel.cr.usgs.gov", 11235);

		// If product consists primarily of binary data, set this option `true`
		// to accelerate distribution.
		sender.setBinaryFormat(false);

		// If product consists primarily of text data, set this option `true`
		// to accelerate distribution.
		sender.setEnableDeflate(true);

		// ^^ Note ^^ Typically do not set both of the above options to `true` as
		//            binary content doesn't compress efficiently but adds
		//            processing overhead.
		
		sender.sendProduct(product);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
