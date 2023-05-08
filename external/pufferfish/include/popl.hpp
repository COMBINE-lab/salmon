/***
    This file is part of popl (program options parser lib)
    Copyright (C) 2015-2016 Johannes Pohl

    This software may be modified and distributed under the terms
    of the MIT license.  See the LICENSE file for details.
***/

#ifndef POPL_H
#define POPL_H

#define NOMINMAX

#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string.h>
#include <vector>

namespace popl {

#define POPL_VERSION "0.3.0"

enum                 // permitted values for its `has_arg' field...
{ no_argument = 0,   // option never takes an argument
  required_argument, // option always requires an argument
  optional_argument  // option may take an argument
};

class Option {
  friend class OptionParser;

public:
  Option(const std::string& shortOption, const std::string& longOption,
         const std::string& description);

  char getShortOption() const;
  std::string getLongOption() const;
  std::string getDescription() const;
  unsigned int count() const;
  bool isSet() const;

protected:
  virtual void parse(const std::string& whatOption, const char* value) = 0;
  virtual void updateReference();
  virtual std::string optionToString() const;
  virtual std::vector<std::string> descriptionToString(size_t width = 40) const;
  virtual int hasArg() const = 0;

  std::string shortOption_;
  std::string longOption_;
  std::string description_;
  unsigned int count_;
};

template <class T> class Value : public Option {
public:
  Value(const std::string& shortOption, const std::string& longOption,
        const std::string& description);
  Value(const std::string& shortOption, const std::string& longOption,
        const std::string& description, const T& defaultVal,
        T* assignTo = NULL);

  Value<T>& assignTo(T* var);
  Value<T>& setDefault(const T& value);
  T getValue(size_t idx = 0) const;

protected:
  virtual void parse(const std::string& whatOption, const char* value);
  virtual std::string optionToString() const;
  virtual int hasArg() const;
  virtual void addValue(const T& value);
  virtual void updateReference();
  T* assignTo_;
  std::vector<T> values_;
  T default_;
  bool hasDefault_;
};

template <class T> class Implicit : public Value<T> {
public:
  Implicit(const std::string& shortOption, const std::string& longOption,
           const std::string& description, const T& implicitVal);
  Implicit(const std::string& shortOption, const std::string& longOption,
           const std::string& description, const T& implicitVal,
           T* assignTo = NULL);

  Value<T>& assignTo(T* var);

protected:
  virtual void parse(const std::string& whatOption, const char* value);
  virtual std::string optionToString() const;
  virtual int hasArg() const;
  Value<T>& setDefault(const T& value);
};

class Switch : public Value<bool> {
public:
  Switch(const std::string& shortOption, const std::string& longOption,
         const std::string& description);
  Switch(const std::string& shortOption, const std::string& longOption,
         const std::string& description, bool* assignTo);

protected:
  virtual void parse(const std::string& whatOption, const char* value);
  virtual std::string optionToString() const;
  virtual int hasArg() const;
  Switch& setDefault(const bool& value);
};

class OptionParser {
public:
  OptionParser(const std::string& description = "");
  virtual ~OptionParser();
  OptionParser& add(Option& option);
  void parse(int argc, char** argv);
  std::string help() const;
  const std::vector<Option*>& options() const;
  const std::vector<std::string>& nonOptionArgs() const;
  const std::vector<std::string>& unknownOptions() const;

protected:
  std::vector<Option*> options_;
  std::string description_;
  std::vector<std::string> nonOptionArgs_;
  std::vector<std::string> unknownOptions_;

  Option* getLongOpt(const std::string& opt) const;
  Option* getShortOpt(char opt) const;
};

/// Option implementation /////////////////////////////////

Option::Option(const std::string& shortOption, const std::string& longOption,
               const std::string& description)
    : shortOption_(shortOption), longOption_(longOption),
      description_(description), count_(0) {
  if (shortOption.size() > 1)
    throw std::invalid_argument("length of short option must be <= 1: '" +
                                shortOption + "'");

  if (shortOption.empty() && longOption.empty())
    throw std::invalid_argument("short and long option are empty");
}

void Option::updateReference() {}

char Option::getShortOption() const {
  if (!shortOption_.empty())
    return shortOption_[0];
  return 0;
}

std::string Option::getLongOption() const { return longOption_; }

std::string Option::getDescription() const { return description_; }

unsigned int Option::count() const { return count_; }

bool Option::isSet() const { return (count() > 0); }

std::string Option::optionToString() const {
  std::stringstream line;
  if (getShortOption() != 0) {
    line << "  -" << getShortOption();
    if (!getLongOption().empty())
      line << ", ";
  } else
    line << "  ";

  if (!getLongOption().empty())
    line << "--" << getLongOption();

  return line.str();
}

std::vector<std::string> Option::descriptionToString(size_t width) const {
  std::vector<std::string> lines;
  std::stringstream description(getDescription());
  std::string line;
  while (std::getline(description, line, '\n'))
    lines.push_back(line);

  return lines;
}

/// Value implementation /////////////////////////////////

template <class T>
Value<T>::Value(const std::string& shortOption, const std::string& longOption,
                const std::string& description)
    : Option(shortOption, longOption, description), assignTo_(NULL),
      hasDefault_(false) {}

template <class T>
Value<T>::Value(const std::string& shortOption, const std::string& longOption,
                const std::string& description, const T& defaultVal,
                T* assignTo)
    : Option(shortOption, longOption, description), assignTo_(assignTo),
      default_(defaultVal), hasDefault_(true) {
  updateReference();
}

template <class T> Value<T>& Value<T>::assignTo(T* var) {
  assignTo_ = var;
  return *this;
}

template <class T> Value<T>& Value<T>::setDefault(const T& value) {
  default_ = value;
  hasDefault_ = true;
  return *this;
}

template <class T> void Value<T>::updateReference() {
  if (assignTo_ != NULL) {
    if (isSet() || hasDefault_)
      *assignTo_ = getValue();
  }
}

template <class T> void Value<T>::addValue(const T& value) {
  values_.push_back(value);
  ++count_;
  updateReference();
}

template <class T> T Value<T>::getValue(size_t idx) const {
  if (!isSet()) {
    if (hasDefault_)
      return default_;
    else {
      std::stringstream optionStr;
      if (getShortOption() != 0)
        optionStr << "-" << getShortOption();
      else
        optionStr << "--" << getLongOption();

      throw std::out_of_range("option not set: \"" + optionStr.str() + "\"");
    }
  }

  if (idx >= count_) {
    std::stringstream optionStr;
    optionStr << "index out of range (" << idx << ") for \"";
    if (getShortOption() != 0)
      optionStr << "-" << getShortOption();
    else
      optionStr << "--" << getLongOption();
    optionStr << "\"";
    throw std::out_of_range(optionStr.str());
  }

  return values_[idx];
}

template <class T> int Value<T>::hasArg() const { return required_argument; }

template <>
void Value<std::string>::parse(const std::string& whatOption,
                               const char* value) {
  if (strlen(value) == 0)
    throw std::invalid_argument("missing argument for " + whatOption);

  addValue(value);
}

template <class T>
void Value<T>::parse(const std::string& whatOption, const char* value) {
  T parsedValue;
  std::string strValue;
  if (value != NULL)
    strValue = value;

  std::istringstream is(strValue);
  int valuesRead = 0;
  while (is.good()) {
    if (is.peek() != EOF)
      is >> parsedValue;
    else
      break;

    valuesRead++;
  }

  if (is.fail())
    throw std::invalid_argument("invalid argument for " + whatOption + ": '" +
                                strValue + "'");

  if (valuesRead > 1)
    throw std::invalid_argument("too many arguments for " + whatOption + ": '" +
                                strValue + "'");

  if (strValue.empty())
    throw std::invalid_argument("missing argument for " + whatOption);

  addValue(parsedValue);
}

template <class T> std::string Value<T>::optionToString() const {
  std::stringstream ss;
  ss << Option::optionToString() << " arg";
  if (hasDefault_) {
    std::stringstream defaultStr;
    defaultStr << default_;
    if (!defaultStr.str().empty())
      ss << " (=" << default_ << ")";
  }
  return ss.str();
}

/// Implicit implementation /////////////////////////////////

template <class T>
Implicit<T>::Implicit(const std::string& shortOption,
                      const std::string& longOption,
                      const std::string& description, const T& implicitVal)
    : Value<T>(shortOption, longOption, description, implicitVal) {}

template <class T>
Implicit<T>::Implicit(const std::string& shortOption,
                      const std::string& longOption,
                      const std::string& description, const T& implicitVal,
                      T* assignTo)
    : Value<T>(shortOption, longOption, description, implicitVal, assignTo) {}

template <class T> int Implicit<T>::hasArg() const { return optional_argument; }

template <class T>
void Implicit<T>::parse(const std::string& whatOption, const char* value) {
  if ((value != NULL) && (strlen(value) > 0))
    Value<T>::parse(whatOption, value);
  else
    this->addValue(this->default_);
}

template <class T> std::string Implicit<T>::optionToString() const {
  std::stringstream ss;
  ss << Option::optionToString() << " [=arg(=" << this->default_ << ")]";
  return ss.str();
}

/// Switch implementation /////////////////////////////////

Switch::Switch(const std::string& shortOption, const std::string& longOption,
               const std::string& description)
    : Value<bool>(shortOption, longOption, description, false) {}

Switch::Switch(const std::string& shortOption, const std::string& longOption,
               const std::string& description, bool* assignTo)
    : Value<bool>(shortOption, longOption, description, false, assignTo) {}

void Switch::parse(const std::string& whatOption, const char* value) {
  addValue(true);
}

int Switch::hasArg() const { return no_argument; }

std::string Switch::optionToString() const { return Option::optionToString(); }

/// OptionParser implementation /////////////////////////////////

OptionParser::OptionParser(const std::string& description)
    : description_(description) {}

OptionParser::~OptionParser() {}

OptionParser& OptionParser::add(Option& option) {
  for (size_t n = 0; n < options_.size(); ++n) {
    if ((option.getShortOption() != 0) &&
        (option.getShortOption() == options_[n]->getShortOption()))
      throw std::invalid_argument("dublicate short option '-" +
                                  std::string(1, option.getShortOption()) +
                                  "'");
    if (!option.getLongOption().empty() &&
        (option.getLongOption() == (options_[n]->getLongOption())))
      throw std::invalid_argument("dublicate long option '--" +
                                  option.getLongOption() + "'");
  }
  options_.push_back(&option);
  return *this;
}

const std::vector<Option*>& OptionParser::options() const { return options_; }

const std::vector<std::string>& OptionParser::nonOptionArgs() const {
  return nonOptionArgs_;
}

const std::vector<std::string>& OptionParser::unknownOptions() const {
  return unknownOptions_;
}

std::string OptionParser::help() const {
  std::stringstream s;
  if (!description_.empty())
    s << description_ << ":\n";

  size_t optionRightMargin(20);
  const size_t maxDescriptionLeftMargin(40);
  //	const size_t descriptionRightMargin(80);

  for (size_t opt = 0; opt < options_.size(); ++opt)
    optionRightMargin =
        std::max(optionRightMargin, options_[opt]->optionToString().size() + 2);
  optionRightMargin = std::min(maxDescriptionLeftMargin - 2, optionRightMargin);

  for (size_t opt = 0; opt < options_.size(); ++opt) {
    std::string optionStr = options_[opt]->optionToString();
    if (optionStr.size() < optionRightMargin)
      optionStr.resize(optionRightMargin, ' ');
    else
      optionStr += "\n" + std::string(optionRightMargin, ' ');
    s << optionStr;

    std::vector<std::string> lines = options_[opt]->descriptionToString(20);
    std::string empty(optionRightMargin, ' ');
    for (size_t n = 0; n < lines.size(); ++n) {
      if (n > 0)
        s << "\n" << empty;
      s << lines[n];
    }
    s << "\n";
  }

  return s.str();
}

Option* OptionParser::getLongOpt(const std::string& opt) const {
  for (size_t n = 0; n < options_.size(); ++n) {
    if (options_[n]->getLongOption() == opt)
      return options_[n];
  }
  return NULL;
}

Option* OptionParser::getShortOpt(char opt) const {
  for (size_t n = 0; n < options_.size(); ++n)
    if (options_[n]->getShortOption() == opt)
      return options_[n];
  return NULL;
}

void OptionParser::parse(int argc, char** argv) {
  unknownOptions_.clear();
  nonOptionArgs_.clear();
  for (int n = 1; n < argc; ++n) {
    const std::string arg(argv[n]);
    if (arg == "--") {
      /// from here on only non opt args
      for (int m = n + 1; m < argc; ++m)
        nonOptionArgs_.push_back(argv[m]);

      break;
    } else if (arg.find("--") == 0) {
      /// long option arg
      std::string opt = arg.substr(2);
      std::string optarg;
      size_t equalIdx = opt.find('=');
      if (equalIdx != std::string::npos) {
        optarg = opt.substr(equalIdx + 1);
        opt.resize(equalIdx);
      }

      Option* option = NULL;
      if ((option = getLongOpt(opt)) != NULL) {
        if (option->hasArg() == no_argument) {
          if (!optarg.empty())
            option = NULL;
        } else if (option->hasArg() == required_argument) {
          if (optarg.empty() && n < argc - 1)
            optarg = argv[++n];
        }
      }

      if (option != NULL)
        option->parse(opt, optarg.c_str());
      else
        unknownOptions_.push_back(arg);
    } else if (arg.find("-") == 0) {
      /// short option arg
      std::string opt = arg.substr(1);
      bool unknown = false;
      for (size_t m = 0; m < opt.size(); ++m) {
        char c = opt[m];
        Option* option = NULL;
        std::string optarg;

        if ((option = getShortOpt(c)) != NULL) {
          if (option->hasArg() == required_argument) {
            /// use the rest of the current argument as optarg
            optarg = opt.substr(m + 1);
            /// or the next arg
            if (optarg.empty() && n < argc - 1)
              optarg = argv[++n];
            m = opt.size();
          } else if (option->hasArg() == optional_argument) {
            /// use the rest of the current argument as optarg
            optarg = opt.substr(m + 1);
            m = opt.size();
          }
        }

        if (option != NULL)
          option->parse(std::string(1, c), optarg.c_str());
        else
          unknown = true;
      }
      if (unknown)
        unknownOptions_.push_back(arg);
    } else {
      nonOptionArgs_.push_back(arg);
    }
  }
}

std::ostream& operator<<(std::ostream& out, const OptionParser& op) {
  out << op.help();
  return out;
}

} // namespace popl

#endif
