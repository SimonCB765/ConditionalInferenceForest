/**
 * 
 */
package tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Simon Bull
 *
 */
public class NonTerminalNumericRegressionNode extends Node
{

	double splitValue;
	Node[] children = new Node[2];
	String covariable;
	Map<String, Double> allVariableTestStats = new HashMap<String, Double>();
	Map<String, Double> allVariablePValues = new HashMap<String, Double>();

	public NonTerminalNumericRegressionNode(String loadString)
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
		this.splitValue = Double.parseDouble(split[5]);
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
	}

	public NonTerminalNumericRegressionNode(int nodeDepth, String covariable, double splitValue,
			Node leftChild, Node rightChild, Map<String, double[]> testStatsAndPVals, int numberOfObservationsInNode)
	{
		this.nodeDepth = nodeDepth;
		this.pValue = 1 + testStatsAndPVals.get(covariable)[2];
		this.splitValue = splitValue;
		this.children[0] = leftChild;
		this.children[1] = rightChild;
		this.covariable = covariable;
		for (String s : testStatsAndPVals.keySet())
		{
			allVariableTestStats.put(s, testStatsAndPVals.get(s)[0]);
			allVariablePValues.put(s, 1 + testStatsAndPVals.get(s)[2]);
		}
		this.numberOfObservationsInNode = numberOfObservationsInNode;
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
			if (!returnValue.contains(this.splitValue))
			{
				returnValue.add(this.splitValue);
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
		System.out.format("Covariable : %s, split value : %f, criterion value %f\n", covariable, splitValue, pValue);
		this.children[0].display();
		this.children[1].display();
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
		double valueOfObservedCovar = (double) currentObservation.get(this.covariable);
		if (valueOfObservedCovar <= this.splitValue)
		{
			return this.children[0].predict(currentObservation);
		}
		else
		{
			return this.children[1].predict(currentObservation);
		}
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
		returnValue += "\t" + Double.toString(this.splitValue) + "\t" + this.covariable + "\t";
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
		return returnValue;
	}

}
