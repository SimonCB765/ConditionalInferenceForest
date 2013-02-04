/**
 * 
 */
package tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Simon Bull
 *
 */
public class NonTerminalCategoricalClassificationNode extends NodeClassification
{

	/**
	 * An array of the same length as the number of categories that there are. The ith entry in splitValue will indicate
	 * the index into children for category i+1 of the covariable. For example:
	 * splitValue[2] == 4, then children.get(4) will be the child node for category 3 of the covariable. 
	 */
	double[] splitValue;
	List<Node> children;
	String covariable;
	Map<String, Double> allVariableTestStats = new HashMap<String, Double>();
	Map<String, Double> allVariablePValues = new HashMap<String, Double>();

	public NonTerminalCategoricalClassificationNode(String loadString)
	{
		loadString = loadString.replaceAll("\n", "");
		String split[] = loadString.split("\t");
		this.nodeDepth = Integer.parseInt(split[0]);
		this.numberOfObservationsInNode = Integer.parseInt(split[1]);
		this.pValue = Double.parseDouble(split[2]);
		for (String s : split[3].split(","))
		{
			String sSplit[] = s.split(";");
			this.pValues.put(sSplit[0], Double.parseDouble(sSplit[1]));
		}
		for (String s : split[4].split(","))
		{
			String sSplit[] = s.split(";");
			this.testStatistics.put(sSplit[0], Double.parseDouble(sSplit[1]));
		}
		String[] splitValueSplit = split[5].split(",");
		this.splitValue = new double[splitValueSplit.length];
		for (int i = 0; i < this.splitValue.length; i++)
		{
			this.splitValue[i] = Double.parseDouble(splitValueSplit[i]);
		}
		this.covariable = split[6];
		for (String s : split[7].split(","))
		{
			String sSplit[] = s.split(";");
			this.allVariableTestStats.put(sSplit[0], Double.parseDouble(sSplit[1]));
		}
		for (String s : split[8].split(","))
		{
			String sSplit[] = s.split(";");
			this.allVariablePValues.put(sSplit[0], Double.parseDouble(sSplit[1]));
		}
		for (String s : split[9].split(","))
		{
			String sSplit[] = s.split(";");
			this.classCountsInNode.put(sSplit[0], Integer.parseInt(sSplit[1]));
		}
	}

	public NonTerminalCategoricalClassificationNode(int nodeDepth, String covariable, double[] splitValue,
			List<Node> children, Map<String, Integer> classCountsInNode, Map<String, double[]> testStatsAndPVals)
	{
		this.nodeDepth = nodeDepth;
		this.pValue = 1 + testStatsAndPVals.get(covariable)[2];
		this.splitValue = splitValue;
		this.children = children;
		this.classCountsInNode = classCountsInNode;
		this.covariable = covariable;
		for (String s : testStatsAndPVals.keySet())
		{
			allVariableTestStats.put(s, testStatsAndPVals.get(s)[0]);
			allVariablePValues.put(s, 1 + testStatsAndPVals.get(s)[2]);
		}
		numberOfObservationsInNode = 0;
		for (String s : this.classCountsInNode.keySet())
		{
			numberOfObservationsInNode += this.classCountsInNode.get(s);
		}
		pValues = new HashMap<String, Double>();
		testStatistics = new HashMap<String, Double>();
		for (String s : testStatsAndPVals.keySet())
		{
			pValues.put(s, 1 + testStatsAndPVals.get(s)[2]);
			testStatistics.put(s, 1 + testStatsAndPVals.get(s)[0]);
		}
	}

	List<Double> cutpoints(String covariable)
	{
		List<Double> returnValue = new ArrayList<Double>();
		for (Node c : this.children)
		{
			returnValue.addAll(c.cutpoints(covariable));
		}
		if (covariable.equals(this.covariable))
		{
			for (double d : this.splitValue)
			{
				returnValue.add(d);
			}
		}
		return returnValue;
	}

	void display()
	{
		for (int i = 0; i < nodeDepth; i++)
		{
			System.out.print("|  ");
		}
		System.out.format("Covariable : %s, splits : ", covariable);
		System.out.println(Arrays.toString(this.splitValue));
		for (Node i : children)
		{
			i.display();
		}
	}

	Set<String> getChildren()
	{
		Set<String> returnValue = new HashSet<String>();
		for (Node c : children)
		{
			returnValue.addAll(c.getChildren());
		}
		returnValue.add(this.covariable);
		return returnValue;
	}

	Object predict(Map<String, Object> currentObservation)
	{
		int valueOfObservedCovar = (int) currentObservation.get(this.covariable);
		double splitIndex = this.splitValue[valueOfObservedCovar - 1];
		Node child = this.children.get((int) splitIndex);
		return child.predict(currentObservation);
	}

	String save()
	{
		String returnValue = "";
		returnValue += Integer.toString(this.nodeDepth) + "/t" + Integer.toString(this.numberOfObservationsInNode) + "\t"
				+ Double.toString(pValue) + "\t";
		for (String s : this.pValues.keySet())
		{
			returnValue += s + ";" + Double.toString(this.pValues.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		for (String s : this.testStatistics.keySet())
		{
			returnValue += s + ";" + Double.toString(this.testStatistics.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		for (Double d : this.splitValue)
		{
			returnValue += Double.toString(d) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		returnValue += this.covariable + "\t";
		for (String s : this.allVariableTestStats.keySet())
		{
			returnValue += s + ";" + Double.toString(this.allVariableTestStats.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		for (String s : this.allVariablePValues.keySet())
		{
			returnValue += s + ";" + Double.toString(this.allVariablePValues.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		returnValue += "\t";
		for (String s : this.classCountsInNode.keySet())
		{
			returnValue += s + ";" + Integer.toString(this.classCountsInNode.get(s)) + ",";
		}
		returnValue = returnValue.substring(0, returnValue.length() - 1);  // Chop off the last ','.
		return returnValue;
	}

}
