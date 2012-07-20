<?php

class RestReader
{
    private $jsObject;

    function __construct($url)
    {
        $file = file_get_contents($url);
        if (!$file) {
            throw new Exception("Forum posts are not loading.");
        } else {
            $this->jsObject = json_decode($file);
            if (!$this->jsObject) {
                throw new Exception("Can't read forum posts.");
            }
        }
    }

    /**
     * Used for usort an array of posts, high to low, based on Views
     */
    private function threadRComp($a, $b)
    {
        // Views is our ranking factor
        $x = $a->CountViews;
        $y = $b->CountViews;

        if ($x == $y) {
            return 0;
        }
        return ($x < $y) ? 1 : -1;
    }

    /**
     * Checks to see of a string contains a particular substring
     * @param $substring the substring to match
     * @param $string the string to search
     * @return true if $substring is found in $string, false otherwise
     */
    private function contains($substring, $string)
    {
        $pos = strpos($string, $substring);
        if ($pos === false) {
            // string needle NOT found in haystack
            return false;
        } else {
            // string needle found in haystack
            return true;
        }
    }

    /**
     * Takes a JSON object and filters all the threads according to a single thread
     * Returns a compatible JSON object
     */
    public function processBy($field, $value)
    {

        foreach ($this->jsObject->Announcements as $post) {
            if ($this->contains($value, $post->$field)) {
                $announcements[] = $post;
            }
        }

        $this->jsObject->Announcements = $announcements;

        foreach ($this->jsObject->Discussions as $post) {
            if ($this->contains($value, $post->$field)) {
                $discussions[] = $post;
            }
        }

        $this->jsObject->Discussions = $discussions;

        return $this;
    }

    private static function printThreads($list, $num, $emptyMsg = "No posts found!")
    {
        echo '<ul class="nav nav-pills nav-stacked">';
        if (count($list) == 0) {
            echo "<li>$emptyMsg</li>";
        } else {
            for ($i = 0; $i < $num && $i < count($list); $i++) {
                $post = $list[$i];
                echo "<li><a href='{$post->Url}'>{$post->Name}</a></li>";
            }
        }
        echo '</ul>';
    }

    /**
     * Prints out an HTML list of the Announcments (by officials) in the given JSON object
     */
    public function printAnnouncements($num = 3, $emptyMsg = "No announcements found")
    {
        echo "<aside class='well'><h2>Official Posts</h2>";
        usort($this->jsObject->Announcements, "threadRComp");
        $this->printThreads($this->jsObject->Announcements, $num, $emptyMsg);
        echo "</aside>";
    }

    /**
     * Prints out an HTML list of the discussions (by users) in the given JSON object
     */
    public function printDiscussions($num = 5, $emptyMsg = "No discussions found")
    {
        echo "<aside class='well'><h2>Popular Discussions</h2>";
        usort($this->jsObject->Discussions, "threadRComp");
        $this->printThreads($this->jsObject->Discussions, $num, $emptyMsg);
        echo "</aside>";
    }


}
