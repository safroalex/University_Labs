package strateg;


import startHero.Point;

public class heStands implements start_stop_point {
    @Override
    public Point newPos(Point start, Point stop) {
        return start;
    }
}
