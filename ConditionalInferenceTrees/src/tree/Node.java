/**
 * 
 */
package tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Simon Bull
 *
 */
public abstract class Node
{

	int nodeDepth;
	int numberOfObservationsInNode;
	double pValue;
	Map<String, Double> pValues;
	Map<String, Double> testStatistics;

	List<Double> cutpoints(String covariable)
	{
		return new ArrayList<Double>();
	}

	void display()
	{
	}

	Set<String> getChildren()
	{
		return new HashSet<String>();
	}

	void load(String loadString)
	{
	}

	Object predict(Map<String, Object> currentObservation)
	{
		return null;
	}

	String save()
	{
		return null;
	}

}
