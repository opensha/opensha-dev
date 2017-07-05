package scratch.aftershockStatistics.cmu;


import javax.jms.Message;

/**
 * Created by clark on 11/11/2016.
 */
public interface MessageListener {
    void onMessage(Message message);
}
