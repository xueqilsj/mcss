/*
 * LogFile.h
 * utility for output log
 *  Created on: Jan 30, 2009
 */

#ifndef LOGFILE_H_
#define LOGFILE_H_

#include <fstream>

#define LOG_STREAM(_log_) (std::fstream &)_log_
#define LOG_INSTANCE(_log_)  _log_.operator std::fstream &()
#define LOG_ERROR(_log_, _what_) LOG_STREAM(_log_) << "!!" << __FILE__ << '(' << __LINE__ << "):\t" << _what_ << std::endl;
#define LOG_ERROR2(_log_, _hint_ , _what_) LOG_STREAM(_log_) << "!!" << __FILE__ << '(' << __LINE__ << "):\t" << _hint_<<_what_ << std::endl;
#define LOG_WARNING(_log_, _what_) LOG_STREAM(_log_) << "%%" << __FILE__ << '(' << __LINE__ << "):\t" << _what_ << std::endl;
#define LOG_MESSAGE(_log_, _what_) LOG_STREAM(_log_) << "**" << __FILE__ << '(' << __LINE__ << "):\t" << _what_ << std::endl;
#define LOG_MESSAGE2(_log_, _hint_ , _what_) LOG_STREAM(_log_) << "**" << __FILE__ << '(' << __LINE__ << "):\t"<< _hint_<<_what_ << std::endl;
#define LOG_INFO(_log_, _name_, _var_, _comment_ ) _log_ << "##" << _name_ << '\t' \
		<< _var_ << '\t'<< '!' << _comment_ << std::endl;

class LogFile {
private:
	///forbid copy!
	LogFile(const LogFile & other);
	LogFile & operator =(const LogFile & other);

protected:
	std::fstream _fStream;
	std::string _SLeadingComment1, _SLeadingComment2; // Reader has get and Writer has set method
public:
	explicit LogFile(const std::string& filename, std::ios::openmode mode =
			std::ios_base::out) {
		_fStream.open(filename.c_str(), mode);
		_SLeadingComment1 = "#";
		_SLeadingComment2 = "!";
	}
	LogFile() {
		_SLeadingComment1 = "#";
		_SLeadingComment2 = "!";
	}
	inline bool Open(const std::string& filename, std::ios::openmode mode =
			std::ios_base::out) {
		if (_fStream.is_open())
			return false;
		_fStream.open(filename.c_str(), mode);
		return _fStream.is_open();
	}
	inline void Close() {
		if (_fStream.is_open())
			_fStream.close();
	}
	inline void SetComment(const std::string &com1, const std::string & com2) {
		_SLeadingComment1 = com1;
		_SLeadingComment2 = com2;
	}
	inline void PutEndl(const std::string &com = "") {
		_fStream << std::endl;
	}
	inline void EatOneLine() {
		char szBuf_at[81];
		_fStream.getline(szBuf_at, 80);
	}
	inline void PutCommentEndl(const std::string &com = "") {
		_fStream << _SLeadingComment2 << com << std::endl;
	}
	inline void PutComment1() {
		_fStream << _SLeadingComment1;
	}
	inline void PutComment2() {
		_fStream << _SLeadingComment2;
	}
	template<class Tvar>
	inline void PutCol(const Tvar & var) {
		_fStream << var << '\t';
	}
	template<class Tvar>
	inline void PutCol(const Tvar & var, const Tvar & var2) {
		_fStream << var << '\t' << var2 << '\t';
	}
	template<class Tvar>
	inline void PutCol(const Tvar & var, const Tvar & var2, const Tvar & var3) {
		_fStream << var << '\t' << var2 << '\t' << var3 << 't';
	}
	template<class Tvar>
	inline void PutInfo(const std::string &name, Tvar & var,
			const std::string & comment) {
		_fStream << _SLeadingComment1 << name << '\t' << var << '\t'
				<< _SLeadingComment2 << comment << std::endl;
	}
	template<class OutputIterator>
	inline void PutRangeInfo(const std::string &name, OutputIterator first,
			OutputIterator last, const std::string & comment) {
		_fStream << _SLeadingComment1 << name << '\t';
		while (first != last)
			_fStream << *first++ << '\t';
		_fStream << _SLeadingComment2 << comment << std::endl;
	}
	template<class OutputIterator>
	inline void PutRange(OutputIterator first, OutputIterator last) {
		while (first != last)
			_fStream << *first++ << '\t';
	}
	template<class Tvar>
	inline void GetCol(Tvar & var) {
		_fStream >> var;
		_fStream.get();
	}
	template<class Tvar>
	inline void GetCol(Tvar & var, Tvar & var2) {
		_fStream >> var >> std::ws >> var2;
		_fStream.get();
	}
	template<class Tvar>
	inline void GetCol(Tvar & var, Tvar & var2, Tvar & var3) {
		_fStream >> var >> std::ws >> var2 >> std::ws >> var3;
		_fStream.get();
	}
	template<class Tvar>
	inline bool GetInfo(const std::string &name, Tvar & var) {
		char __szBuf[81];
		_fStream.getline(__szBuf, 80, '\t');
		if (name.compare(__szBuf + _SLeadingComment1.length()) != 0)
			return false;
		_fStream >> var;
		_fStream.getline(__szBuf, 80);
		return true;
	}
	template<class InputIterator>
	inline bool GetRangeInfo(const std::string &name, InputIterator first,
			InputIterator last) {
		char szBuf_at[81];
		_fStream.getline(szBuf_at, 80, '\t');
		if (name.compare(szBuf_at + _SLeadingComment1.length()) != 0)
			return false;
		while (first != last) {
			_fStream >> *(first++);//< not std::ios_base::skipws
			_fStream.get();//<skip '\t'
		}
		_fStream.getline(szBuf_at, 80);
		return true;
	}

	template<class InputIterator>
	inline void GetRange(InputIterator first, InputIterator last) {
		while (first != last) {
			_fStream >> *(first++);//< not std::ios_base::skipws
			_fStream.get();//<skip '\t'
		}
	}
	inline bool IsOpen() {
		return _fStream.is_open();
	}
	///(std::fstream &) for <<,>>
	inline operator std::fstream &() {
		return _fStream;
	}
	~LogFile() {
		if (_fStream.is_open())
			_fStream.close();

	}

};

#endif /* LOGFILE_H_ */
