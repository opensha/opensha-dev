package scratch.aftershockStatistics.cmu;

import javax.jms.*;
import javax.jms.MessageConsumer;

import org.apache.activemq.ActiveMQConnectionFactory;

/**
 * Created by clark on 11/11/2016.
 */
public class Jobs {
    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration();

    	System.out.println("Jobs Service Starting...");
        final ConnectionFactory connectionFactory = new ActiveMQConnectionFactory("tcp://" + config.getActivemq_host() + ":" + config.getActivemq_port());
        final Connection connection = connectionFactory.createConnection(config.getActivemq_user(), config.getActivemq_password());

        Runtime.getRuntime().addShutdownHook(new Thread(){
            @Override
            public void run(){
                try {

                	System.out.println("Shutdown Hook Activated...");
                    if(connection != null){
                        connection.close();
                    }
                } catch (JMSException e) {
                    e.printStackTrace();
                }
            }
        });

        try {
            final Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);

            /** Jobs.RJ Consumer */
            Queue queue = session.createQueue("Jobs.RJ");
            MessageConsumer consumer = session.createConsumer(queue);
            consumer.setMessageListener(new JobsListener("Jobs Queue Consumer"));
            connection.start();

        	System.out.println("Jobs Service Listening to the Queue Jobs.RJ...");
            while(true){
                Thread.sleep(Integer.MAX_VALUE);
            }
        } finally {
        	System.out.println("Jobs Service Exiting...");
            if (connection != null) {
                connection.close();
            }
        }
    }
}
