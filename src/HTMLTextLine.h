/*
 * Header file for HTMLTextLine
 * Copyright (C) 2013 Lu Wang <coolwanglu@gmail.com>
 */
#ifndef HTMLTEXTLINE_H__
#define HTMLTEXTLINE_H__

#include <ostream>
#include <vector>
#include <list>

#include <CharTypes.h>

#include "Param.h"
#include "StateManager.h"
#include "HTMLState.h"

namespace pdf2htmlEX {

/*
 * Store and optimize a line of text in HTML
 *
 * contains a series of 
 *  - Text
 *  - Shift
 *  - State change
 */
class HTMLTextLine
{
public:
    //Added by Tyler Clemens.
    struct LetterState {
        Unicode* letter;  //the characters
        int length;       //the length of the characters
        double x, y;      //character position
        double dx, dy;    //character width and height
        int index;        //the index of the character on its line
        double fs;        //Font Size
        double dts;       //Draw Text Scale
    };

    //Added by Tyler Clemens
    struct WordState {
        std::list<LetterState*>::iterator first_letter, last_letter;
        double avgls;                       // average letter space >= 0. round negative values up to 0.
        double x, y;                        // the coordinates of the word
        double mcu_cs;                      // the most commonly used space between characters. Exclude negatives and 0
        double als;                         // the average letter space

        void print(std::ostream &out);
    };

    HTMLTextLine (const HTMLLineState & line_state, const Param & param, AllStateManager & all_manager);
    ~HTMLTextLine();

    struct State : public HTMLTextState {
        State();
        // before output
        void begin(std::ostream & out, const State * prev_state);
        // after output
        void end(std::ostream & out) const;
        // calculate the hash code
        void hash(void);
        // calculate the difference between another State
        int diff(const State & s) const;

        //Added by Tyler Clemens.
        void append_letter_state(std::list<LetterState*>::iterator end){
            std::list<LetterState*>::iterator last = end;
            last--;
            if(begining)
                first_letter = last;
            last_letter = last;
            begining = false;
        }
        void set_mcu_cs();
        std::list<WordState>::iterator detect_spaces_and_split(std::list<WordState>::iterator beginWord, std::list<WordState>::iterator word);

        enum {
            FONT_ID,
            FONT_SIZE_ID,
            FILL_COLOR_ID,
            STROKE_COLOR_ID,
            LETTER_SPACE_ID,
            WORD_SPACE_ID,
            HASH_ID_COUNT,

            VERTICAL_ALIGN_ID = HASH_ID_COUNT,
            ID_COUNT
        };

        static long long umask_by_id(int id);

        long long ids[ID_COUNT];

        size_t start_idx; // index of the first Text using this state
        // for optimzation
        long long hash_value;
        long long hash_umask; // some states may not be actually used
        bool need_close;

        static const char * const css_class_names []; // class names for each id

        //Added by Tyler Clemens. A vector of letter states
        //std::list<LetterState*> letters;
        std::list<LetterState*>::iterator first_letter, last_letter;
        //Added by Tyler Clemens. A vector of words
        std::list<WordState> words;
        double x, y;
        double dts;
        std::map<double, int> cses;
        double mcu_cs;                  // most commonly used letter space
        bool begining;
    
    };

    struct Offset {
        Offset(size_t size_idx, double width)
            :start_idx(size_idx),width(width)
        { }
        size_t start_idx; // should put this Offset right before text[start_idx];
        double width;
    };

    void append_unicodes(const Unicode * u, int l);
    void append_offset(double width);
    void append_state(const HTMLTextState & text_state);

    // Added by Tyler Clemens. A method for appending to LetterPositions
    void append_letter_state(Unicode *letter, int uLen, double x, double y, double dx, double dy, int index, double fs, double dts); 
    void make_words(std::list<State>::iterator cur_state);
    
    void dump_text(std::ostream & out);

    bool text_empty(void) const { return text.empty(); }
    void clear(void);

    void clip(const HTMLClipState &);

    /*
     * Optimize and calculate necessary values
     */
    void prepare(void);
    std::list<State> states;
    std::list<LetterState*> letters;
private:
    void optimize(void);

    //Added by Tyler Clemens. A quick and dirty way to output words when given a string of characters
    void outputWords(std::ostream & out, const Unicode * u, int uLen, std::list<State>::iterator state1, int cur_text_idx);
    void calculateWordPos(std::vector<LetterState>::iterator, double &x, double &y);
    bool checkForSpace(std::vector<LetterState>::iterator letters);

    const Param & param;
    AllStateManager & all_manager;

    HTMLLineState line_state;
    double ascent, descent;
    double clip_x1, clip_y1;

    std::vector<Offset> offsets;
    std::vector<Unicode> text;

    //Added by Tyler Clemens.
    double all_spaces;
    double avg_char_space;
    std::list<State>::iterator cur_state_itr;
    double x, y; // the x and y position of the state on the line
};

} // namespace pdf2htmlEX
#endif //HTMLTEXTLINE_H__
