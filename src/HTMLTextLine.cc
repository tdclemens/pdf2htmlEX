/*
 * HTMLTextLine.cc
 *
 * Generate and optimized HTML for one line
 *
 * Copyright (C) 2012,2013 Lu Wang <coolwanglu@gmail.com>
 */

#include <cmath>
#include <algorithm>

#include "HTMLTextLine.h"

#include "util/encoding.h"
#include "util/css_const.h"

namespace pdf2htmlEX {

using std::min;
using std::max;
using std::vector;
using std::ostream;
using std::cerr;
using std::endl;
using std::find;
using std::abs;

HTMLTextLine::HTMLTextLine (const HTMLLineState & line_state, const Param & param, AllStateManager & all_manager) 
    :param(param)
    ,all_manager(all_manager) 
    ,line_state(line_state)
    ,clip_x1(0)
    ,clip_y1(0)
{ 
}

HTMLTextLine::~HTMLTextLine(){
    for(std::list<LetterState*>::iterator itr = letters.begin(); itr != letters.end(); itr++){
        delete[] (*itr)->letter;
        delete *itr;
    }
}

void HTMLTextLine::append_unicodes(const Unicode * u, int l)
{
    text.insert(text.end(), u, u+l);
}

void HTMLTextLine::append_offset(double width)
{
    /*
     * If the last offset is very thin, we can ignore it and directly use it
     * But this should not happen often, and we will also filter near-zero offsets when outputting them
     * So don't check it
     */
    if((!offsets.empty()) && (offsets.back().start_idx == text.size()))
        offsets.back().width += width;
    else
        offsets.emplace_back(text.size(), width);
}

void HTMLTextLine::append_state(const HTMLTextState & text_state)
{
    if(states.empty() || (states.back().start_idx != text.size()))
    {
        states.emplace_back();
        states.back().start_idx = text.size();
        states.back().hash_umask = 0;
    }

    HTMLTextState & last_state = states.back();
    last_state = text_state;
    //apply font scale
    last_state.font_size *= last_state.font_info->font_size_scale;
}


void HTMLTextLine::append_letter_state(Unicode *letter, int uLen, double x, double y, double dx, double dy, int index, double fs, double dts) 
{
    letters.push_back(new LetterState());
    //letters.back()->letter.insert(letters.back()->letter.end(), letter, letter + uLen);
    letters.back()->letter = new Unicode[uLen];
    for(int i = 0; i < uLen; i++)
        (letters.back()->letter)[i] = letter[i];
    letters.back()->length = uLen;
    letters.back()->x = x;
    letters.back()->dx = dx;
    letters.back()->dy = dy;
    letters.back()->index = index;
    letters.back()->fs = fs;
    letters.back()->dts = dts;
}


void HTMLTextLine::dump_text(ostream & out){

    //std::ostream &out = std::cout;

    //set the values of x for all letters
    {
      std::list<LetterState*>::iterator cur_letter, next_letter;
      double prevx = (*letters.begin())->x;
//      if(prevx < 0)
//          printf("prevx: %f\n", prevx);
      //letters.begin()->x = letters.begin()->x * letters.begin()->dts;
      (*(letters.begin()))->x = 0;
      cur_letter = letters.begin();
      next_letter = letters.begin();
      next_letter++;
      for(; next_letter != letters.end(); cur_letter++, next_letter++){
          double x  = ((*next_letter)->x - prevx) * (*cur_letter)->dts;
          prevx = (*next_letter)->x;
          (*next_letter)->x = (*cur_letter)->x + x;
//          if(next_letter->x < 0)
//              printf("%f\n", next_letter->x);
      }
      /*
              std::list<LetterState>::iterator cur_letter, next_letter, unscaled_cur, unscaled_next;
              cur_letter = next_letter = letters.begin();
          if(letters.size() > 0)
              cur_letter->x = cur_letter->x * cur_letter->dts;
          if(letters.size() > 1){
              std::list<LetterState> unscaledLetters;
              std::copy(letters.begin(), letters.end(), unscaledLetters.begin());
              printf("after copy\n");
              unscaled_cur = unscaled_next = unscaledLetters.begin();
              next_letter++;
              unscaled_next++;
              for(; next_letter != letters.end(); cur_letter++, next_letter++){
                  next_letter->x = (unscaled_next->x - unscaled_cur->x) * unscaled_cur->dts;
                  if(next_letter->x < 0)
                      printf("%f\n",next_letter->x);
              }
         }
         */
      //set the first state on the line to have the right x coordinate
      states.front().x = (*states.front().first_letter)->x;
      //loop through each state
      for(std::list<State>::iterator cur_state = states.begin(); cur_state != states.end(); cur_state++){
          //set values for each state
          cur_state->dts = (*(cur_state->first_letter))->dts;
          cur_state->x = (*(cur_state->first_letter))->x;
          //gather the characters into words based on space
          make_words(cur_state);
      }
    }

    //set the most common character space of each state
    for(std::list<State>::iterator cur_state = states.begin(); cur_state != states.end(); cur_state++)
        cur_state->set_mcu_cs();

    //for each state, try to guess where there is a space and insert new words
    std::list<WordState>::iterator cur_word;
    for(std::list<State>::iterator cur_state = states.begin(); cur_state != states.end(); cur_state++){
//            cur_state->detect_spaces_and_split(cur_state->words.begin(),cur_state->words.begin());
    }
    
    // do some of the original stuff
    if(states.empty() || (states.front().start_idx != 0))
    {
        cerr << "Warning: text without a style! Must be a bug in pdf2htmlEX" << endl;
        return;
    }

    // Start Output
    {
        // open <div> for the current text line
        out << "<div class=\"" << CSS::LINE_CN
            << " " << CSS::TRANSFORM_MATRIX_CN << all_manager.transform_matrix.install(line_state.transform_matrix)
            << " " << CSS::LEFT_CN             << all_manager.left.install(line_state.x - clip_x1)
            << " " << CSS::HEIGHT_CN           << all_manager.height.install(ascent)
            << " " << CSS::BOTTOM_CN           << all_manager.bottom.install(line_state.y - clip_y1)
            ;
        // it will be closed by the first state
    }
    State * prev_state = nullptr;
    for(std::list<State>::iterator cur_state = states.begin(); cur_state != states.end(); cur_state++){
        cur_state->begin(out, prev_state);
        for(std::list<WordState>::iterator words = cur_state->words.begin(); words != cur_state->words.end(); words++)
            words->print(out);
        cur_state->end(out);
        prev_state = &(*cur_state);
    }

    out << "</div>";

}

void HTMLTextLine::make_words(std::list<State>::iterator cur_state){
    //loop through each character in the state and make a word if the current character is a space.
    cur_state->words.emplace_back();
    cur_state->words.back().x = (*(cur_state->first_letter))->x;
    cur_state->words.back().als = 0;
    cur_state->words.back().first_letter = cur_state->first_letter;
    std::list<LetterState*>::iterator copy, prev_letter, end_last_letter;;
    prev_letter = cur_state->first_letter;
    end_last_letter = cur_state->last_letter;
    end_last_letter++;
    int ct = 0;
    bool begining = true;
    for(std::list<LetterState*>::iterator cur_letter = cur_state->first_letter; cur_letter != end_last_letter; cur_letter++){
        cur_state->words.back().last_letter = cur_letter;
        if(!begining)
            cur_state->words.back().als += (*cur_letter)->x - (*prev_letter)->x - (*prev_letter)->dx * (*prev_letter)->dts;
        //search the letter string for a space
        bool found = false;
        for(int i = 0; i < (*cur_letter)->length; i++)
            //if(((*cur_letter)->letter)[i] == ' ')
            if(((*cur_letter)->letter)[i] == ' ' || ((*cur_letter)->letter)[i] == ',' || ((*cur_letter)->letter)[i] == '.' || ((*cur_letter)->letter)[i] == ';' || ((*cur_letter)->letter)[i] == ':')
                found = true;
        //if you find a space create a new word
        if(found){
            cur_state->words.back().als = cur_state->words.back().als /( ct - 1);
            ct = 0;
            copy = cur_letter;
            copy++;
            cur_state->words.emplace_back();
            if(copy != end_last_letter){
                cur_state->words.back().x = (*copy)->x;
                cur_state->words.back().first_letter = copy;
            }
            else{
                cur_state->words.back().x = (*cur_letter)->x;
                cur_state->words.back().first_letter = cur_state->words.back().last_letter = cur_letter;
            }
            //prev_letter = cur_state->letters.begin();
            cur_state->words.back().als = 0;
        }
        else{
            ct++;
            prev_letter = cur_letter;
      }
      begining = false;
    }
}

void HTMLTextLine::State::set_mcu_cs(){
    std::list<LetterState*>::iterator cur_letter, next_letter;
    std::list<LetterState*>::iterator endItr;
    for(std::list<WordState>::iterator cur_word = words.begin(); cur_word != words.end(); cur_word++){
        endItr = cur_word->last_letter;;
        //endItr--;
        //TODO: change this algorithm to a more efficient one. This is a brute force method
        if(std::distance(cur_word->first_letter, cur_word->last_letter) > 1)
        for(cur_letter = cur_word->first_letter,next_letter = cur_letter,next_letter++; next_letter != endItr; cur_letter++, next_letter++){
            double cs = (*next_letter)->x - (*cur_letter)->x - ((*cur_letter)->dx * (*cur_letter)->dts);
            if(cs > 0){
                if((cses).find(cs) == cses.end()){
                    cses[cs] = 1;
                }
                else{
                    cses[cs] = cses[cs] + 1;
                }
            }
        }
    }
    mcu_cs = 0;
    for(std::map<double,int>::iterator itr = cses.begin(); itr != cses.end(); itr++){
        if(mcu_cs < itr->first)
            mcu_cs = itr->first;
    }
    //printf("state mcu_cs: %f\n", mcu_cs);

}

std::list<HTMLTextLine::WordState>::iterator HTMLTextLine::State::detect_spaces_and_split(std::list<WordState>::iterator beginWord, std::list<WordState>::iterator word){
    if(word == words.end())
        return word;
    // end case for recursion
    std::ostream &out = std::cout;

    std::list<LetterState*>::iterator cur_letter, next_letter, end_last_letter;
    std::list<WordState>::iterator copy = word;
    copy++;
    WordState newWord;
    bool inserted = false;
    end_last_letter = word->last_letter;
    end_last_letter++;
    if(std::distance(word->first_letter, word->last_letter) > 1){
        for(cur_letter = word->first_letter,next_letter = cur_letter, next_letter++; next_letter != end_last_letter; cur_letter++, next_letter++){
            double cs = (*next_letter)->x - (*cur_letter)->x - ((*cur_letter)->dx * (*cur_letter)->dts);
            if(cs * dts >= single_space_offset()){
                //set the x value for the word
                newWord.x = (*next_letter)->x;
                //copy the rest of the characters over to the new word and delete them
                newWord.first_letter = next_letter;
                for(std::list<LetterState*>::iterator itr = next_letter; itr != end_last_letter; itr++)
                    newWord.last_letter = itr;
                //word->letters.erase(next_letter, word->letters.end());
                //insert the new word before the next word
                words.insert(copy, newWord);
                break;
            }
        }
        copy = word;
        copy++;
        return detect_spaces_and_split(beginWord, copy);
    }
    else
      return word;
}

void HTMLTextLine::clear(void)
{
    states.clear();
    offsets.clear();
    text.clear();
}

void HTMLTextLine::clip(const HTMLClipState & clip_state)
{
    clip_x1 = clip_state.xmin;
    clip_y1 = clip_state.ymin;
}

void HTMLTextLine::prepare(void)
{
    if(param.optimize_text)
        optimize();

    // max_ascent determines the height of the div
    double accum_vertical_align = 0; // accumulated
    ascent = 0;
    descent = 0;
    // note that vertical_align cannot be calculated here
    for(auto iter = states.begin(); iter != states.end(); ++iter)
    {
        auto font_info = iter->font_info;
        iter->ids[State::FONT_ID] = font_info->id;
        iter->ids[State::FONT_SIZE_ID]      = all_manager.font_size.install(iter->font_size);
        iter->ids[State::FILL_COLOR_ID]     = all_manager.fill_color.install(iter->fill_color);
        iter->ids[State::STROKE_COLOR_ID]   = all_manager.stroke_color.install(iter->stroke_color);
        iter->ids[State::LETTER_SPACE_ID]   = all_manager.letter_space.install(iter->letter_space);
        iter->ids[State::WORD_SPACE_ID]     = all_manager.word_space.install(iter->word_space);
        iter->hash();

        accum_vertical_align += iter->vertical_align;
        double cur_ascent = accum_vertical_align + font_info->ascent * iter->font_size;
        if(cur_ascent > ascent)
            ascent = cur_ascent;
        double cur_descent = accum_vertical_align + font_info->descent * iter->font_size;
        if(cur_descent < descent)
            descent = cur_descent;
    }
}


/*
 * Adjust letter space and word space in order to reduce the number of HTML elements
 * May also unmask word space
 */
void HTMLTextLine::optimize()
{
    // remove unuseful states in the end
    while((!states.empty()) && (states.back().start_idx >= text.size()))
        states.pop_back();

    assert(!states.empty());

    const long long word_space_umask = State::umask_by_id(State::WORD_SPACE_ID);

    // for optimization, we need accurate values
    auto & ls_manager = all_manager.letter_space;
    auto & ws_manager = all_manager.word_space;
    
    // statistics of widths
    std::map<double, size_t> width_map;
    // store optimized offsets
    std::vector<Offset> new_offsets;
    new_offsets.reserve(offsets.size());

    auto offset_iter1 = offsets.begin();
    for(auto state_iter2 = states.begin(), state_iter1 = state_iter2++; 
            state_iter1 != states.end(); 
            ++state_iter1, ++state_iter2)
    {
        const size_t text_idx1 = state_iter1->start_idx;
        const size_t text_idx2 = (state_iter2 == states.end()) ? text.size() : state_iter2->start_idx;
        size_t text_count = text_idx2 - text_idx1;

        // there might be some offsets before the first state
        while((offset_iter1 != offsets.end()) 
                && (offset_iter1->start_idx <= text_idx1))
        {
            new_offsets.push_back(*(offset_iter1++));
        }

        // find the last offset covered by the current state
        auto offset_iter2 = offset_iter1;
        for(; (offset_iter2 != offsets.end()) && (offset_iter2->start_idx <= text_idx2); ++offset_iter2) { }

        // There are `offset_count` <span>'s, the target is to reduce this number
        size_t offset_count = offset_iter2 - offset_iter1;
        assert(text_count >= offset_count);
        
        // Optimize letter space
        // how much letter_space is changed
        // will be later used for optimizing word space
        double letter_space_diff = 0; 
        width_map.clear();

        // In some PDF files all letter spaces are implemented as position shifts between each letter
        // try to simplify it with a proper letter space
        if(offset_count > 0)
        {
            // mark the current letter_space
            if(text_count > offset_count)
                width_map.insert(std::make_pair(0, text_count - offset_count));

            for(auto off_iter = offset_iter1; off_iter != offset_iter2; ++off_iter)
            {
                const double target = off_iter->width;
                auto iter = width_map.lower_bound(target-EPS);
                if((iter != width_map.end()) && (abs(iter->first - target) <= EPS))
                {
                    ++ iter->second;
                }
                else
                {
                    width_map.insert(iter, std::make_pair(target, 1));
                }
            }
            
            // TODO snapping the widths may result a better result
            // e.g. for (-0.7 0.6 -0.2 0.3 10 10), 0 is better than 10
            double most_used_width = 0;
            size_t max_count = 0;
            for(auto iter = width_map.begin(); iter != width_map.end(); ++iter)
            {
                if(iter->second > max_count)
                {
                    most_used_width = iter->first;
                    max_count = iter->second;
                }
            }

            // negative letter space may cause problems
            if((max_count <= text_count / 2) || (!is_positive(state_iter1->letter_space + most_used_width)))
            { 
                // the old value is the best
                // just copy old offsets
                new_offsets.insert(new_offsets.end(), offset_iter1, offset_iter2);
            }
            else
            {
                // now we would like to adjust letter space to most_used width
                
                // install new letter space
                const double old_ls = state_iter1->letter_space;
                state_iter1->ids[State::LETTER_SPACE_ID] = ls_manager.install(old_ls + most_used_width, &(state_iter1->letter_space));
                letter_space_diff = old_ls - state_iter1->letter_space;
                // update offsets
                auto off_iter = offset_iter1; 
                // re-count number of offsets
                offset_count = 0;
                for(size_t cur_text_idx = text_idx1; cur_text_idx < text_idx2; ++cur_text_idx)
                {
                    double cur_width = 0;
                    if((off_iter != offset_iter2) && (off_iter->start_idx == cur_text_idx + 1))
                    {
                        cur_width = off_iter->width + letter_space_diff;
                        ++off_iter;
                    }
                    else
                    {
                        cur_width = letter_space_diff ;
                    }
                    if(!equal(cur_width, 0))
                    {
                        new_offsets.emplace_back(cur_text_idx+1, cur_width);
                        ++ offset_count;
                    }
                }
            }
        }

        // Optimize word space
        
        // In some PDF files all spaces are converted into positionig shift
        // We may try to change (some of) them to ' ' by adjusting word_space
        // for now, we cosider only the no-space scenario
        // which also includes the case when param.space_as_offset is set

        // get the text segment covered by current state (*state_iter1)
        const auto text_iter1 = text.begin() + text_idx1;
        const auto text_iter2 = text.begin() + text_idx2;
        if(find(text_iter1, text_iter2, ' ') == text_iter2)
        {
            // if there is not any space, we may change the value of word_space arbitrarily
            // note that we may only change word space, no offset will be affected
            // The actual effect will emerge during flushing, where it could be detected that an offset can be optimized as a single space character
            
            if(offset_count > 0)
            {
                double threshold = (state_iter1->em_size()) * (param.space_threshold);
                // set word_space for the most frequently used offset
                double most_used_width = 0;
                size_t max_count = 0;

                // if offset_count > 0, we must have updated width_map in the previous step
                // find the most frequent width, with new letter space applied
                for(auto iter = width_map.begin(); iter != width_map.end(); ++iter)
                {
                    double fixed_width = iter->first + letter_space_diff; // this is the actual offset in HTML
                    // we don't want to add spaces for tiny gaps, or even negative shifts
                    if((fixed_width >= threshold - EPS) && (iter->second > max_count))
                    {
                        max_count = iter->second;
                        most_used_width = fixed_width;
                    }
                }

                state_iter1->word_space = 0; // clear word_space for single_space_offset
                double new_word_space = most_used_width - state_iter1->single_space_offset();
                state_iter1->ids[State::WORD_SPACE_ID] = ws_manager.install(new_word_space, &(state_iter1->word_space)); // install new word_space
                state_iter1->hash_umask &= (~word_space_umask); // mark that the word_space is not free
            }
            else // there is no offset at all
            {
                state_iter1->hash_umask |= word_space_umask; // we just free word_space
            }
        }
        offset_iter1 = offset_iter2;
    } 
    
    // apply optimization
    std::swap(offsets, new_offsets);
}

// this state will be converted to a child node of the node of prev_state
// dump the difference between previous state
// also clone corresponding states
void HTMLTextLine::State::begin (ostream & out, const State * prev_state)
{
    if(prev_state)
    {
        long long cur_mask = 0xff;
        bool first = true;
        for(int i = 0; i < HASH_ID_COUNT; ++i, cur_mask<<=8)
        {
            if(hash_umask & cur_mask) // we don't care about this ID
            {
                if (prev_state->hash_umask & cur_mask) // if prev_state do not care about it either
                    continue;

                // otherwise
                // we have to inherit it
                ids[i] = prev_state->ids[i]; 
                hash_umask &= (~cur_mask);
                //copy the corresponding value
                //TODO: this is so ugly
                switch(i)
                {
                    case FONT_SIZE_ID:
                        font_size = prev_state->font_size;
                        break;
                    case LETTER_SPACE_ID:
                        letter_space = prev_state->letter_space;
                        break;
                    case WORD_SPACE_ID:
                        word_space = prev_state->word_space;
                        break;
                    default:
                        cerr << "unexpected state mask" << endl;
                        break;
                }
            }

            // now we care about the ID
            
            // if the value from prev_state is the same, we don't need to dump it
            if((!(prev_state->hash_umask & cur_mask)) && (prev_state->ids[i] == ids[i]))
                continue;

            // so we have to dump it
            if(first)
            { 
                out << "<span style='position:absolute;bottom:0px;height:inherit;left:0px;' class=\"";
                first = false;
            }
            else
            {
                out << ' ';
            }

            // out should have hex set
            out << css_class_names[i];
            if (ids[i] == -1)
                out << CSS::INVALID_ID;
            else
                out << ids[i];
        }
        // veritcal align
        if(!equal(vertical_align, 0))
        {
            // so we have to dump it
            if(first)
            { 
                out << "<span style='position:absolute;bottom:0px;height:inherit;left:0px;' class=\"";
                first = false;
            }
            else
            {
                out << ' ';
            }

            // out should have hex set
            out << CSS::VERTICAL_ALIGN_CN;
            auto id = ids[VERTICAL_ALIGN_ID];
            if (id == -1)
                out << CSS::INVALID_ID;
            else
                out << id;
        }

        if(first) // we actually just inherit the whole prev_state
        {
            need_close = false;
        }
        else
        {
            out << "\">";
            need_close = true;
        }
    }
    else
    {
        // prev_state == nullptr
        // which means this is the first state of the line
        // there should be a open pending <div> left there
        // it is not necessary to output vertical align
        long long cur_mask = 0xff;
        for(int i = 0; i < HASH_ID_COUNT; ++i, cur_mask<<=8)
        {
            if(hash_umask & cur_mask) // we don't care about this ID
                continue;

            // now we care about the ID
            out << ' '; 
            // out should have hex set
            out << css_class_names[i];
            if (ids[i] == -1)
                out << CSS::INVALID_ID;
            else
                out << ids[i];
        }

        out << "\">";
        need_close = false;
    }
}

void HTMLTextLine::State::end(ostream & out) const
{
    if(need_close)
        out << "</span>";
}

void HTMLTextLine::State::hash(void)
{
    hash_value = 0;
    for(int i = 0; i < ID_COUNT; ++i)
    {
        hash_value = (hash_value << 8) | (ids[i] & 0xff);
    }
}

int HTMLTextLine::State::diff(const State & s) const
{
    /*
     * A quick check based on hash_value
     * it could be wrong when there are more then 256 classes, 
     * in which case the output may not be optimal, but still 'correct' in terms of HTML
     */
    long long common_mask = ~(hash_umask | s.hash_umask);
    if((hash_value & common_mask) == (s.hash_value & common_mask)) return 0;

    long long cur_mask = 0xff;
    int d = 0;
    for(int i = 0; i < ID_COUNT; ++i)
    {
        if((common_mask & cur_mask) && (ids[i] != s.ids[i]))
            ++ d;
        cur_mask <<= 8;
    }
    return d;
}

long long HTMLTextLine::State::umask_by_id(int id)
{
    return (((long long)0xff) << (8*id));
}

HTMLTextLine::State::State(){
    begining = true;
}

//print the word with its beginning and ending span
void HTMLTextLine::WordState::print(std::ostream &out){
  std::list<LetterState*>::iterator backItr = last_letter;
  last_letter++;
  out << "<span style='position:absolute;bottom:0px;left:"<<x<<"px;letter-spacing:"<<"0"<<"px;line-height:50%;'>";
  for(std::list<LetterState*>::iterator itr = first_letter; itr != last_letter; itr++){
      for(int i = 0; i < (*itr)->length; i++){
          outputUnicodes(out, &(((*itr)->letter)[i]), 1);
          //insert a space at the end of the word if its not there
          if(itr == backItr &&  i == (*itr)->length - 1 && ((*itr)->letter)[i] != ' ')
              out << ' ';
      }
  }
  out << "</span>";
}

// the order should be the same as in the enum
const char * const HTMLTextLine::State::css_class_names [] = {
    CSS::FONT_FAMILY_CN,
    CSS::FONT_SIZE_CN,
    CSS::FILL_COLOR_CN,
    CSS::STROKE_COLOR_CN,
    CSS::LETTER_SPACE_CN,
    CSS::WORD_SPACE_CN,
    CSS::VERTICAL_ALIGN_CN,
};

} //namespace pdf2htmlEX
