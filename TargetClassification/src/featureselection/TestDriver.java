/**
 * 
 */
package featureselection;

import tree.TreeGrowthControl;

/**
 * @author Simon
 *
 */
public class TestDriver
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		TreeGrowthControl ctrl = new TreeGrowthControl();
		ctrl.minCriterion = -Double.MAX_VALUE;
		ctrl.isClassificationUsed = true;
		ctrl.isReplacementUsed = true;
		ctrl.numberOfTreesToGrow = 100;
		int gaRepetitions = 10;
		boolean isXValUsed = false;
		new Controller(args, ctrl, gaRepetitions, isXValUsed);
	}

}
