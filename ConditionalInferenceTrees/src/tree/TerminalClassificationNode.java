/**
 * 
 */
package tree;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Simon Bull
 *
 */
public class TerminalClassificationNode extends NodeClassification
{

	Map<String, Double> classWeightsInNode;

	TerminalClassificationNode(String loadString)
	{
		loadString = loadString.replaceAll("\n", "");
		String split[] = loadString.split("\t");
		this.classCountsInNode = new HashMap<String, Integer>();
		for (String s : split[0].split(","))
		{
			String sSplit[] = s.split(";");
			this.classCountsInNode.put(sSplit[0], Integer.parseInt(sSplit[1]));
		}
		this.nodeDepth = Integer.parseInt(split[1]);
		this.pValue = 0.0;
		this.pValues = new HashMap<String, Double>();
		this.testStatistics = new HashMap<String, Double>();
		for (String s : this.classCountsInNode.keySet())
		{
			this.pValues.put(s, 0.0);
			this.testStatistics.put(s, 0.0);
		}
	}

	TerminalClassificationNode(Map<String, Integer> classCounts, int currentDepth, Map<String, Double> classWeightsInNode)
	{
		this.classCountsInNode = classCounts;
		this.nodeDepth = currentDepth;
		this.classWeightsInNode = classWeightsInNode;
		this.pValue = 0.0;
		pValues = new HashMap<String, Double>();
		testStatistics = new HashMap<String, Double>();
		for (String s : classCounts.keySet())
		{
			pValues.put(s, 0.0);
			testStatistics.put(s, 0.0);
		}
	}

	void display()
	{
		for (int i = 0; i < nodeDepth; i++)
		{
			System.out.print("|  ");
		}
		System.out.println(classCountsInNode.entrySet());
	}

	Object predict(Map<String, Object> currentObservation)
	{
		Map<String, Double> classFractions = new HashMap<String, Double>();
		double numberObservationsInNode = 0.0;
		for (Integer j : classCountsInNode.values())
		{
			numberObservationsInNode += j;
		}
		for (String className : classCountsInNode.keySet())
		{
			classFractions.put(className, classCountsInNode.get(className) / numberObservationsInNode);
		}
		Object returnValue = this.classWeightsInNode;

		return returnValue;
	}

	String save()
	{
		String returnValue = "";
		for (String s : this.classCountsInNode.keySet())
		{
			returnValue += s + ";" + Integer.toBinaryString(this.classCountsInNode.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		returnValue += Integer.toString(this.nodeDepth);
		return returnValue;
	}

}
