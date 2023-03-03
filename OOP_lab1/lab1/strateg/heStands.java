package strateg;


import common.Point;

public class heStands implements MoveStrategy {
    @Override
    public Point getNewPos(Point start, Point stop) {
        return start;
    }
}
